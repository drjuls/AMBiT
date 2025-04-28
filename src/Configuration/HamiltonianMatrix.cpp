#include "Include.h"
#include "HamiltonianMatrix.h"
#include "HartreeFock/Orbital.h"
#include "Projection.h"
#include "Universal/Eigensolver.h"
#include "Universal/MathConstant.h"
#include "Universal/ScalapackMatrix.h"
#include <gsl/gsl_statistics_ulong.h>
#include <memory>
#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif

#ifdef AMBIT_USE_OPENMP
#include<omp.h>
#endif


// Don't bother with davidson method if smaller than this limit
#define SMALL_MATRIX_LIM 200

// Don't bother with davidson method if number of solutions requested is larger than this
// and ScaLAPACK is available.
#define MANY_LEVELS_LIM   50

namespace Ambit
{
HamiltonianMatrix::HamiltonianMatrix(pHFIntegrals hf, pTwoElectronCoulombOperator coulomb, pRelativisticConfigList relconfigs):
    H_two_body(nullptr), H_three_body(nullptr), configs(relconfigs), most_chunk_rows(0)
{
    // Set up Hamiltonian operator
    H_two_body = std::make_shared<TwoBodyHamiltonianOperator>(hf, coulomb);

    // Set up matrix
    N = configs->NumCSFs();
    Nsmall = configs->NumCSFsSmall();

    if(Nsmall != N)
    {
        *logstream << " " << N << "x" << Nsmall << std::endl;
        *outstream << " Number of CSFs = " << N << " x " << Nsmall << std::flush;
    }
    else
    {
        *logstream << " " << N << " " << std::flush;
        *outstream << " Number of CSFs = " << N << std::flush;
    }
}

HamiltonianMatrix::HamiltonianMatrix(pHFIntegrals hf, pTwoElectronCoulombOperator coulomb, pSigma3Calculator sigma3, pConfigListConst leadconfigs, pRelativisticConfigList relconfigs):
    HamiltonianMatrix(hf, coulomb, relconfigs)
{
    // Set up three-body operator
    H_three_body = std::make_shared<ThreeBodyHamiltonianOperator>(hf, coulomb, sigma3);
    leading_configs = leadconfigs;
}

HamiltonianMatrix::~HamiltonianMatrix()
{}

void HamiltonianMatrix::GenerateMatrix(unsigned int configs_per_chunk)
{
    chunks.clear();

    if(N <= SMALL_MATRIX_LIM)
    {
        configs_per_chunk = configs->size();
    }

    // Total number of chunks = ceiling(number of configs/configs_per_chunk)
    unsigned int total_num_chunks = (configs->size() + configs_per_chunk - 1)/configs_per_chunk;

    // Divide up chunks
    auto config_it = configs->begin();
    unsigned int config_index = 0;
    unsigned int csf_start = 0;

    // Loop through the chunks but don't actually construct any yet. We just want to work out the
    // load-balancing among MPI ranks and assign chunks to ranks
    std::vector<size_t> chunks_work_sizes; // Work for each chunk
    for(int chunk_index = 0; chunk_index < total_num_chunks; chunk_index++)
    {
        // Get chunk num_rows and number of configs. Allocates resources for the chunks.
        unsigned int current_num_rows = 0;
        unsigned int current_num_configs = 0;
        size_t current_chunk_work_units = 0;
        while(config_it != configs->end() && current_num_configs < configs_per_chunk)
        {
            current_chunk_work_units += config_it->projection_size()*config_it->projection_size()*config_it->NumCSFs();
            current_num_rows += config_it->NumCSFs();
            current_num_configs++;
            config_it++;
        }

        if(current_num_rows == 0)
            break;

        chunks_work_sizes.push_back(current_chunk_work_units);
    }

    /* Note:
     * Now work out some statistics about the distribution of work among chunks. We want to
     * get the median amount of work per chunk, as well as the Croux-Rousseuw Qn measure of
     * spread (which is more robust when dealing with highly-skewed distributions like this one)
     * to identify chunks which are particularly bad for workload balancing. See this paper for
     * more info on this measure:
     * https://wis.kuleuven.be/stat/robust/papers/publications-1993/rousseeuwcroux-alternativestomedianad-jasa-1993.pdf
     *
     * Also see the GSL manual for information on how it's implemented in GSL:
     * https://www.gnu.org/software/gsl/doc/html/statistics.html#robust-scale-estimates
     */
    double median, Qn;

    // New scope to ensure tmp arrays are deallocated ASAP
    {
      // Temporary workspace arrays for GSL
      std::vector<size_t> gsl_work(3*total_num_chunks);
      std::vector<int> gsl_work_int(5*total_num_chunks);

      // Need to make a deep copy of the data since it must be sorted in ascending order for GSL to
      // calculate the stats
      std::vector<size_t> tmp = chunks_work_sizes;
      std::sort(tmp.begin(), tmp.end());
      // GSL needs a raw pointer to the work size data
      size_t* worksize_pointer = tmp.data();

      median = gsl_stats_ulong_median_from_sorted_data(tmp.data(), 1, total_num_chunks);

      /* Note:
       * The magic factor of 1.566 is necessary here because Qn includes a magic weighting
       * factor based on the assumed distribution of the data. GSL uses the magic constant for a
       * Gaussian distribution, but the workload data is EXTREMELY non-Gaussian due to its skew.
       * It sort of looks exponential if you squint at it, so that's the value I'm using here
       */
      Qn = 1.566*gsl_stats_ulong_Qn_from_sorted_data(tmp.data(), 1, total_num_chunks, gsl_work.data(), gsl_work_int.data());
    }

    /* Note:
     * Now calculate the outlier threshold: any chunk with more than median + 9.0*Qn work units is
     * considered to be a "big chunk". Note that this threshold is somewhat arbitrary, but seems to
     * work okay (think of it as an analogy to the 1.5*IQR rule, but designed for highly-skewed
     * data)
    */
    size_t outlier_threshold = median + 9.0*Qn;

    // Now do another passthrough and actually construct this rank's chunks
    std::vector<size_t> processor_work_sizes(NumProcessors, 0); // Work assigned to each MPI rank
    config_it = configs->begin();
    config_index = 0;
    csf_start = 0;
    int num_big_chunks = 0;

    for(int chunk_index = 0; chunk_index < total_num_chunks; chunk_index++)
    {
        // Get chunk num_rows and number of configs. Allocates resources for the chunks.
        unsigned int current_num_rows = 0;
        unsigned int current_num_configs = 0;
        size_t current_chunk_work_units = 0;
        while(config_it != configs->end() && current_num_configs < configs_per_chunk)
        {
            current_chunk_work_units += config_it->projection_size()*config_it->projection_size()*config_it->NumCSFs();
            current_num_rows += config_it->NumCSFs();
            current_num_configs++;
            config_it++;
        }

        if(current_num_rows == 0)
            break;

        // Assign this chunk to whichever process currently has the least work
        auto min_work_it = std::min_element(processor_work_sizes.begin(), processor_work_sizes.end());
        int assigned_process = std::distance(processor_work_sizes.begin(), min_work_it);
        // Now make the chunk if it's ours
        if(assigned_process == ProcessorRank)
        {
            bool is_big_chunk; 
            if (current_chunk_work_units >= outlier_threshold){
                is_big_chunk = true;
                num_big_chunks++;
            }
            else
            {
                is_big_chunk = false;
            }
            chunks.emplace_back(config_index, config_index+current_num_configs, csf_start,
                                current_num_rows, Nsmall, is_big_chunk);
        }
        // This needs to be outside the conditional so it gets executed by each rank. Every process
        // needs to know how much work has already been assigned to the others
        processor_work_sizes[assigned_process] += current_chunk_work_units;

        config_index += current_num_configs;
        csf_start += current_num_rows;
        most_chunk_rows = mmax(most_chunk_rows, current_num_rows);
    }
    // Print some diagnostics about the workload balancing
    // What is this process's chunk workload?
    *logstream << "This process has " << processor_work_sizes[ProcessorRank] << " work units and " << num_big_chunks << " big chunks" << std::endl;

    // How big is the workload imbalance?
    auto min_work_it = std::min_element(processor_work_sizes.begin(), processor_work_sizes.end());
    auto max_work_it = std::max_element(processor_work_sizes.begin(), processor_work_sizes.end());
    double imbalance = 100.0*((double) (*max_work_it) - (double) (*min_work_it))/((double) (*min_work_it));

    *logstream << "Minimum workload: " << *min_work_it << std::endl;
    *logstream << "Maximum workload: " << *max_work_it << std::endl;
    *logstream << "The relative workload imbalance across MPI processes is " << imbalance << "%" << std::endl;
    // Loop through my chunks
    RelativisticConfigList::const_iterator configsubsetend_it = configs->small_end();
    unsigned int configsubsetend = configs->small_size();

#ifdef AMBIT_USE_OPENMP
    #pragma omp parallel private(config_it)
    {
    #pragma omp single nowait
    {
    #pragma omp taskloop
#endif
    for(size_t chunk_index = 0; chunk_index < chunks.size(); chunk_index++)
    {
        auto& current_chunk = chunks[chunk_index];
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& M = current_chunk.chunk;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& D = current_chunk.diagonal;

        // Loop through configs for this chunk
        config_it = (*configs)[current_chunk.config_indices.first];
        for(unsigned int config_index = current_chunk.config_indices.first; config_index < current_chunk.config_indices.second; config_index++)
        {
            bool leading_config_i = H_three_body && std::binary_search(leading_configs->first.begin(), leading_configs->first.end(), NonRelConfiguration(*config_it));

            // Loop through the rest of the configs
            auto config_jt = configs->begin();
            RelativisticConfigList::const_iterator config_jend;
            if(config_index < configsubsetend)
            {   config_jend = config_it;
                ++config_jend;
            }
            else
                config_jend = configsubsetend_it;

#ifdef AMBIT_USE_OPENMP
            #pragma omp task untied \
                             default(none) \
                             shared(M) \
                             firstprivate(leading_config_i, config_it, config_jt, config_jend, \
                                          current_chunk)
            {
#endif
            while(config_jt != config_jend)
            {
                bool leading_config_j = H_three_body && std::binary_search(leading_configs->first.begin(), leading_configs->first.end(), NonRelConfiguration(*config_jt));

                int config_diff_num = config_it->GetConfigDifferencesCount(*config_jt);
                bool do_three_body = (leading_config_i || leading_config_j) && (config_diff_num <= 3);

                // Check that the number of differences is small enough
                if(do_three_body || (config_diff_num <= 2))
                {
                    // Loop through projections
                    auto proj_it = config_it.projection_begin();
                    while(proj_it != config_it.projection_end())
                    {
                        RelativisticConfiguration::const_projection_iterator proj_jt;
                        if(config_jt == config_it)
                            proj_jt = proj_it;
                        else
                            proj_jt = config_jt.projection_begin();

                        while(proj_jt != config_jt.projection_end())
                        {
                            double operatorH;
                            if(do_three_body)
                            {
                                operatorH = H_three_body->GetMatrixElement(*proj_it, *proj_jt);
                            }
                            else
                            {
                                operatorH = H_two_body->GetMatrixElement(*proj_it, *proj_jt);
                            }
                            if(fabs(operatorH) > 1.e-15)
                            {
                                for(auto coeff_i = proj_it.CSF_begin(); coeff_i != proj_it.CSF_end(); coeff_i++)
                                {
                                    RelativisticConfigList::const_CSF_iterator start_j = proj_jt.CSF_begin();

                                    if(proj_it == proj_jt)
                                        start_j = coeff_i;

                                    for(auto coeff_j = start_j; coeff_j != proj_jt.CSF_end(); coeff_j++)
                                    {
                                        // See notes for an explanation
                                        int i = coeff_i.index();
                                        int j = coeff_j.index();

                                        if(i > j)
                                            M(i - current_chunk.start_row, j) += operatorH * (*coeff_i) * (*coeff_j);
                                        else if(i < j)
                                            M(j - current_chunk.start_row, i) += operatorH * (*coeff_i) * (*coeff_j);
                                        else if(proj_it == proj_jt)
                                            M(i - current_chunk.start_row, j) += operatorH * (*coeff_i) * (*coeff_j);
                                        else
                                            M(i - current_chunk.start_row, j) += 2. * operatorH * (*coeff_i) * (*coeff_j);
                                    }
                                }
                            }
                            proj_jt++;
                        }
                        proj_it++;
                    }
                }
                config_jt++;
            } // while (config_jt)
#ifdef AMBIT_USE_OPENMP
            } // OMP task
#endif

            // Diagonal
            if(config_index >= configs->small_size())
            {
#ifdef AMBIT_USE_OPENMP
                #pragma omp task untied \
                                 default(none) \
                                 shared(D) \
                                 firstprivate(current_chunk, config_it)
                {
#endif
                int diag_offset = current_chunk.start_row + current_chunk.num_rows - current_chunk.diagonal.rows();

                // Loop through projections
                auto proj_it = config_it.projection_begin();
                while(proj_it != config_it.projection_end())
                {
                    RelativisticConfiguration::const_projection_iterator proj_jt = proj_it;

                    while(proj_jt != config_it.projection_end())
                    {
                        double operatorH = H_two_body->GetMatrixElement(*proj_it, *proj_jt);

                        if(fabs(operatorH) > 1.e-15)
                        {
                            for(auto coeff_i = proj_it.CSF_begin(); coeff_i != proj_it.CSF_end(); coeff_i++)
                            {
                                RelativisticConfigList::const_CSF_iterator start_j = proj_jt.CSF_begin();

                                if(proj_it == proj_jt)
                                    start_j = coeff_i;

                                for(auto coeff_j = start_j; coeff_j != proj_jt.CSF_end(); coeff_j++)
                                {
                                    // See notes for an explanation
                                    int i = coeff_i.index();
                                    int j = coeff_j.index();

                                    if(i > j)
                                        D(i - diag_offset, j - diag_offset) += operatorH * (*coeff_i) * (*coeff_j);
                                    else if(i < j)
                                        D(j - diag_offset, i - diag_offset) += operatorH * (*coeff_i) * (*coeff_j);
                                    else if(proj_it == proj_jt)
                                        D(i - diag_offset, j - diag_offset) += operatorH * (*coeff_i) * (*coeff_j);
                                    else
                                        D(i - diag_offset, j - diag_offset) += 2. * operatorH * (*coeff_i) * (*coeff_j);
                                }
                            }
                        }
                        proj_jt++;
                    }
                    proj_it++;
                } // while (proj_it)
#ifdef AMBIT_USE_OPENMP
                } // OMP task
#endif
            } // Conditional for diagonal elements
            config_it++;
        } // Configs in chunk
    } // Chunks loop
#ifdef AMBIT_USE_OPENMP
    } // OMP single
    }  // OMP Parallel
#endif

    for(auto& matrix_section: chunks)
        matrix_section.Symmetrize();
}

LevelVector HamiltonianMatrix::SolveMatrix(pHamiltonianID hID, unsigned int num_solutions)
{
    LevelVector levelvec(hID);
    levelvec.configs = configs;

    unsigned int NumSolutions = mmin(num_solutions, N);

    if(NumSolutions == 0)
    {
        *outstream << "\nNo solutions" << std::endl;
    }
    else
    {   if(N <= SMALL_MATRIX_LIM && NumProcessors == 1 && Nsmall == N)
        {
            *outstream << "; Finding solutions using Eigen..." << std::endl;
            levelvec.levels.reserve(NumSolutions);

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(chunks.front().chunk);
            const Eigen::VectorXd& E = es.eigenvalues();
            const Eigen::MatrixXd& V = es.eigenvectors();

            for(unsigned int i = 0; i < NumSolutions; i++)
            {
                levelvec.levels.push_back(std::make_shared<Level>(E(i), V.col(i).data(), hID, N));
            }
        }
        else if(NumSolutions > MANY_LEVELS_LIM && Nsmall == N)
        {
            *outstream << "; Attempting to reallocate matrix and find solutions using Eigen..." << std::endl;
            levelvec.levels.reserve(NumSolutions);

            RowMajorMatrix M = RowMajorMatrix::Zero(N, N);
            for(auto& chunk: chunks)
            {
                M.block(chunk.start_row, 0, chunk.chunk.rows(), chunk.chunk.cols()) = chunk.chunk;
            }

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
            const Eigen::VectorXd& E = es.eigenvalues();
            const Eigen::MatrixXd& V = es.eigenvectors();

            for(unsigned int i = 0; i < NumSolutions; i++)
            {
                levelvec.levels.push_back(std::make_shared<Level>(E(i), V.col(i).data(), hID, N));
            }
        }
        else
        {   *outstream << "; Finding solutions using Davidson..." << std::endl;
            levelvec.levels.reserve(NumSolutions);

            double* V = new double[NumSolutions * N];
            double* E = new double[NumSolutions];

            Eigensolver solver;
            #ifdef AMBIT_USE_MPI
                solver.MPISolveLargeSymmetric(this, E, V, N, NumSolutions);
            #else
                solver.SolveLargeSymmetric(this, E, V, N, NumSolutions);
            #endif

            for(unsigned int i = 0; i < NumSolutions; i++)
            {
                levelvec.levels.push_back(std::make_shared<Level>(E[i], (V + N * i), hID, N));
            }

            delete[] E;
            delete[] V;
        }
    }

    return levelvec;
}

#ifdef AMBIT_USE_SCALAPACK
LevelVector HamiltonianMatrix::SolveMatrixScalapack(pHamiltonianID hID, unsigned int num_solutions, bool use_energy_limit, double energy_limit)
{
    LevelVector levelvec;
    levelvec.hID = hID;
    levelvec.configs = configs;

    unsigned int NumSolutions = mmin(num_solutions, N);

    if(NumSolutions == 0)
    {
        *outstream << "\nNo solutions" << std::endl;
    }
    else
    {   *outstream << "; Finding solutions using ScaLAPACK ..." << std::endl;

        // Write temporary matrix file, clear current Hamiltonian to make space,
        // then read in to ScalapackMatrix
        char* jobid = getenv("PBS_JOBID");
        std::string filename = "temp";
        if(jobid)
            filename += jobid;
        filename += ".matrix";

        Write(filename);
        Clear();

        ScalapackMatrix SM(N);
        SM.ReadLowerTriangle(filename);

        // Diagonalise
        double* E = new double[N];  // All eigenvalues
        SM.Diagonalise(E);

        // Cut off num_solutions
        if(use_energy_limit)
            for(int i = 0; i < NumSolutions; i++)
            {
                if(E[i] > energy_limit)
                {
                    NumSolutions = i;
                    break;
                }
            };

        // Get levels. Using a larger buffer is generally better, so
        // choose something around 1M * 8 bytes (small enough to be "in the noise")
        levelvec.levels.reserve(NumSolutions);
        unsigned int column_begin = 0;
        unsigned int num_columns_per_step = 1000000/N;
        double* V = new double[N * num_columns_per_step];  // Eigenvectors

        while(column_begin < NumSolutions)
        {
            unsigned int column_end = mmin(column_begin + num_columns_per_step, NumSolutions);
            SM.GetColumns(column_begin, column_end, V);

            double* pV = V;
            while(column_begin < column_end)
            {   levelvec.levels.push_back(std::make_shared<Level>(E[column_begin], pV, hID, N));
                column_begin++;
                pV += N;
            }
        }

        delete[] E;
        delete[] V;
    }

    return levelvec;
}
#endif

std::ostream& operator<<(std::ostream& stream, const HamiltonianMatrix& matrix)
{
    for(auto& matrix_section: matrix.chunks)
    {
        // Each row separately
        for(unsigned int row = 0; row < matrix_section.num_rows; row++)
        {
            int cols = mmin(matrix_section.start_row + row + 1, matrix.Nsmall);

            // Lower triangular matrix part of row
            stream << matrix_section.chunk.block(row, 0, 1, cols) << " ";

            // Trailing zeros
            stream << Eigen::VectorXd::Zero(matrix.Nsmall - cols).transpose() << "\n";
        }
    }
    stream.flush();

    return stream;
}

bool HamiltonianMatrix::Read(const std::string& filename)
{
    FILE* fp = file_err_handler->fopen(filename.c_str(), "rb");
    if(!fp)
        return false;

    // Read size of matrix
    unsigned int matrix_N;
    file_err_handler->fread(&matrix_N, sizeof(unsigned int), 1, fp);

    if(matrix_N != N)
    {   *errstream << "HamiltonianMatrix::Read: matrix has incorrect dimension (number of CSFs)." << std::endl;
        file_err_handler->fclose(fp);
        return false;
    }

    auto chunk_it = chunks.begin();
    double buf[N];
    Eigen::Map<Eigen::RowVectorXd> buf_map(buf, N);

    for(unsigned int matrix_row = 0; matrix_row < N && chunk_it != chunks.end(); matrix_row++)
    {
        file_err_handler->fread(buf, sizeof(double), matrix_row+1, fp);

        // Find chunk for this row, if it exists in this process
        while(matrix_row >= chunk_it->start_row + chunk_it->num_rows && chunk_it != chunks.end())
            chunk_it++;
        if(chunk_it == chunks.end())
            break;
        if(matrix_row < chunk_it->start_row)
            continue;

        // Change Map array using the "placement new" syntax.
        // This does not invoke the memory allocator, because the syntax specifies the location for storing the result.
        new (&buf_map) Eigen::Map<Eigen::RowVectorXd>(buf,matrix_row+1);

        // Copy row into chunk
        unsigned int row_length = mmin(matrix_row+1, Nsmall);
        chunk_it->chunk.block(matrix_row-chunk_it->start_row, 0, 1, row_length) = buf_map.head(row_length);

        // Copy row section into diagonal
        if(Nsmall <= matrix_row)
        {
            // Row in diagonal
            unsigned int diag_row = chunk_it->diagonal.rows() - (chunk_it->num_rows + chunk_it->start_row - matrix_row);
            // Number of elements to copy = diag_row+1;
            chunk_it->diagonal.block(diag_row, 0, 1, diag_row+1) = buf_map.tail(diag_row+1);
        }
    }

    file_err_handler->fclose(fp);

    for(auto& chunk: chunks)
        chunk.Symmetrize();

    return true;
}

void HamiltonianMatrix::Write(const std::string& filename) const
{
    FILE* fp;

    auto chunk_it = chunks.begin();

    // Send rows to root node, which writes them sequentially.
    if(ProcessorRank == 0)
    {
        // Write size of matrix.
        fp = file_err_handler->fopen(filename.c_str(), "wb");
        file_err_handler->fwrite(&N, sizeof(unsigned int), 1, fp);
        const double* pbuf;
        const double* pdiag;
        std::vector<double> zeros(N-Nsmall, 0.);

    #ifdef AMBIT_USE_MPI
        double buf[Nsmall * most_chunk_rows];
        double diagbuf[most_chunk_rows * most_chunk_rows];
    #endif

        int row = 0;
        while(row < N)
        {
            int num_rows = 0;
            int diag_rows = 0;

            // My chunk!
            if(row == chunk_it->start_row)
            {
                num_rows = chunk_it->num_rows;
                pbuf = chunk_it->chunk.data();
                diag_rows = chunk_it->diagonal.rows();
                pdiag = chunk_it->diagonal.data();
                chunk_it++;
            }
        #ifdef AMBIT_USE_MPI
            else
            {   // Broadcast row number
                MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);

                // Receive chunk
                MPI_Status status;
                MPI_Recv(&buf, Nsmall*most_chunk_rows, MPI_DOUBLE, MPI_ANY_SOURCE, row, MPI_COMM_WORLD, &status);

                // Get number of rows in chunk
                int data_count;
                MPI_Get_count(&status, MPI_DOUBLE, &data_count);
                if(data_count >= int(Nsmall) * (int(Nsmall) - row))
                    num_rows = data_count/Nsmall;
                else
                    num_rows = (-row + sqrt(row * row + 4 * data_count))/2;

                if(num_rows * mmin(row + num_rows, Nsmall) != data_count)
                {   *errstream << "HamiltonianMatrix::Write: received incorrect chunk size." << std::endl;
                    exit(1);
                }

                pbuf = buf;

                // Receive diagonal
                if(row + num_rows > Nsmall)
                {
                    MPI_Recv(&diagbuf, most_chunk_rows*most_chunk_rows, MPI_DOUBLE, MPI_ANY_SOURCE, row+1, MPI_COMM_WORLD, &status);

                    // Check diagonal size
                    diag_rows = mmin(num_rows, row + num_rows - Nsmall);
                    MPI_Get_count(&status, MPI_DOUBLE, &data_count);
                    if(diag_rows * diag_rows != data_count)
                    {   *errstream << "HamiltonianMatrix::Write: received incorrect diagonal chunk size." << std::endl;
                        exit(1);
                    }

                    pdiag = diagbuf;
                }
            }
        #endif

            int gap = mmax(row - int(Nsmall), 0);

            for(int i = 0; i < num_rows; i++)
            {
                // Write chunk
                file_err_handler->fwrite(pbuf, sizeof(double), mmin(row + i + 1, Nsmall), fp);
                pbuf += mmin(row + num_rows, Nsmall);   // Move to next row

                // Write diagonal
                if(diag_rows && row + i >= Nsmall)
                {
                    if(gap)
                        file_err_handler->fwrite(zeros.data(), sizeof(double), gap, fp);
                    file_err_handler->fwrite(pdiag, sizeof(double), row + i + 1 - Nsmall - gap, fp);
                    pdiag += diag_rows;
                }
            }

            row += num_rows;
        }

    #ifdef AMBIT_USE_MPI
        // Send finished signal
        MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);
    #endif

        file_err_handler->fclose(fp);
    }
    else
    {
    #ifdef AMBIT_USE_MPI
        int row = 0;
        while(row < N)
        {
            // Received row number
            MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // If it is our row, send chunk
            if(chunk_it != chunks.end() && row == chunk_it->start_row)
            {
                MPI_Send(chunk_it->chunk.data(), chunk_it->chunk.size(), MPI_DOUBLE, 0, row, MPI_COMM_WORLD);

                // Send diagonal if it exists
                if(chunk_it->diagonal.size())
                    MPI_Send(chunk_it->diagonal.data(), chunk_it->diagonal.size(), MPI_DOUBLE, 0, row+1, MPI_COMM_WORLD);

                chunk_it++;
            }
        }
    #endif
    }

#ifdef AMBIT_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

double HamiltonianMatrix::PollMatrix(double epsilon) const
{
    unsigned int i, j, count = 0;
    double value;

    // Iterate over chunks
    for(const auto& it: chunks)
    {
        for(i = 0; i < it.num_rows; i++)
            for(j = i; j < it.chunk.cols(); j++)
                if(fabs(it.chunk(i, j)) > epsilon)
                    count++;
    }

    value = double(count)/double(N * (N+1)/2);
    return value;
}

void HamiltonianMatrix::MatrixMultiply(int m, double* b, double* c) const
{
    Eigen::Map<Eigen::MatrixXd> b_mapped(b, N, m);
    Eigen::Map<Eigen::MatrixXd> c_mapped(c, N, m);
    c_mapped = Eigen::MatrixXd::Zero(N, m);

    // Multiply each chunk
    for(const auto& matrix_section: chunks)
    {
        unsigned int start = matrix_section.start_row;
        unsigned int cols = matrix_section.chunk.cols();

        // Lower triangular part
        c_mapped.middleRows(start, matrix_section.num_rows)
            += matrix_section.chunk * b_mapped.topRows(cols);

        // Upper triangular part
        if(start > 0)
        {   unsigned int upper1_rows = mmin(start, Nsmall);
            c_mapped.topRows(upper1_rows)
                += matrix_section.chunk.leftCols(upper1_rows).transpose() * b_mapped.middleRows(start, matrix_section.num_rows);
        }

        // Extra upper part
        if(start < Nsmall && Nsmall < (start + matrix_section.num_rows))
        {
            unsigned int upper2_cols = start + matrix_section.num_rows - Nsmall;

            c_mapped.middleRows(start, Nsmall - start)
                += matrix_section.chunk.block(Nsmall - start, start, upper2_cols, Nsmall - start).transpose()
                    * b_mapped.middleRows(Nsmall, upper2_cols);
        }

        // Diagonal part
        if(matrix_section.diagonal.rows())
        {
            unsigned int diag_rows  = matrix_section.diagonal.rows();
            unsigned int diag_start = matrix_section.start_row + matrix_section.num_rows - diag_rows;
            c_mapped.middleRows(diag_start, diag_rows)
                += matrix_section.diagonal * b_mapped.middleRows(diag_start, diag_rows);
        }
    }
}

void HamiltonianMatrix::GetDiagonal(double* diag) const
{
    // Get diagonal from each chunk
    Eigen::Map<Eigen::VectorXd> diag_mapped(diag, N, 1);
    diag_mapped.noalias() = Eigen::VectorXd::Zero(N);

    for(const auto& matrix_section: chunks)
    {
        if(matrix_section.start_row < Nsmall)
        {
            unsigned int length = mmin(matrix_section.num_rows, Nsmall - matrix_section.start_row);
            diag_mapped.segment(matrix_section.start_row, length).noalias()
                = matrix_section.chunk.rightCols(length).diagonal();
        }

        if(Nsmall < matrix_section.start_row + matrix_section.num_rows)
        {
            unsigned int start = matrix_section.start_row + matrix_section.num_rows - matrix_section.diagonal.rows();
            diag_mapped.segment(start, matrix_section.diagonal.rows())
                = matrix_section.diagonal.diagonal();
        }
    }
}
}
