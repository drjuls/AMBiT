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

// Don't bother with splitting matrix if smaller than this limit
#define SMALL_MATRIX_LIM 200

// Don't bother with davidson method if number of solutions requested is larger than this
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
    {   // Place whole matrix in a single chunk
        if(ProcessorRank == 0)
            chunks.emplace_back(0, configs->size(), 0, N, Nsmall);
    }
    else
    {   // Estimate workloads per configuration
        // We just want to work out the load-balancing among MPI ranks and assign chunks to ranks
        std::vector<size_t> config_work;    // Work for each configuration
        config_work.reserve(configs->size());
        auto config_it = configs->begin();
        while (config_it != configs->end())
        {
            config_work.push_back(config_it->projection_size() * config_it->projection_size() * config_it->NumCSFs());
            config_it++;
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
            std::vector<size_t> gsl_work(3 * config_work.size());
            std::vector<int> gsl_work_int(5 * config_work.size());

            // Need to make a deep copy of the data since it must be sorted in ascending order for GSL to
            // calculate the stats. We reverse the ordering to begin with since the config list is probably already
            // sorted by workload
            std::vector<size_t> tmp(config_work.rbegin(), config_work.rend());
            std::sort(tmp.begin(), tmp.end());
            // GSL needs a raw pointer to the work size data
            size_t *worksize_pointer = tmp.data();

            median = gsl_stats_ulong_median_from_sorted_data(tmp.data(), 1, config_work.size());

            /* Note:
             * The magic factor of 1.566 is necessary here because Qn includes a magic weighting
             * factor based on the assumed distribution of the data. GSL uses the magic constant for a
             * Gaussian distribution, but the workload data is EXTREMELY non-Gaussian due to its skew.
             * It sort of looks exponential if you squint at it, so that's the value I'm using here
             */
            Qn = 1.566 * gsl_stats_ulong_Qn_from_sorted_data(tmp.data(), 1, config_work.size(), gsl_work.data(),
                                                             gsl_work_int.data());
        }

        /* Note:
         * Now calculate the outlier threshold: any chunk with more than median + 9.0*Qn work units is
         * considered to be a "big chunk". Note that this threshold is somewhat arbitrary, but seems to
         * work okay (think of it as an analogy to the 1.5*IQR rule, but designed for highly-skewed
         * data)
        */
        size_t outlier_threshold = median + 9.0 * Qn;

        // Now do another passthrough and actually construct this rank's chunks
        std::vector<size_t> processor_work_sizes(NumProcessors, 0); // Work assigned to each MPI rank
        unsigned int config_index = 0;
        unsigned int csf_start = 0;
        int num_big_chunks = 0;

        // Start with big chunks, one config at a time
        config_it = configs->begin();
        while (config_it != configs->end())
        {
            if (config_work[config_index] >= outlier_threshold)
            {
                // Assign this chunk to whichever process currently has the least work
                auto min_work_it = std::min_element(processor_work_sizes.begin(), processor_work_sizes.end());
                int assigned_process = std::distance(processor_work_sizes.begin(), min_work_it);
                // Now make the chunk if it's ours
                if (assigned_process == ProcessorRank)
                {
                    num_big_chunks++;
                    chunks.emplace_back(config_index, config_index + 1, csf_start, config_it->NumCSFs(), Nsmall);
                }

                // This needs to be outside the conditional so it gets executed by each rank. Every process
                // needs to know how much work has already been assigned to the others
                processor_work_sizes[assigned_process] += config_work[config_index];
                most_chunk_rows = mmax(most_chunk_rows, config_it->NumCSFs());
            }
            csf_start += config_it->NumCSFs();
            config_index++;
            config_it++;
        }

        // Do the rest of the chunks according to configs_per_chunk suggestion
        config_it = configs->begin();
        config_index = 0;
        csf_start = 0;
        while (config_it != configs->end())
        {
            if (config_work[config_index] >= outlier_threshold)
            {
                csf_start += config_it->NumCSFs();
                config_index++;
                config_it++;
                continue;
            }

            unsigned int current_num_rows = 0;
            unsigned int current_num_configs = 0;
            size_t current_chunk_work_units = 0;
            while (config_it != configs->end() && current_num_configs < configs_per_chunk
                   && config_work[config_index + current_num_configs] < outlier_threshold)
            {
                current_chunk_work_units += config_work[config_index + current_num_configs];
                current_num_rows += config_it->NumCSFs();
                current_num_configs++;
                config_it++;
            }

            if (current_num_rows > 0)
            {   // Assign this chunk to whichever process currently has the least work
                auto min_work_it = std::min_element(processor_work_sizes.begin(), processor_work_sizes.end());
                int assigned_process = std::distance(processor_work_sizes.begin(), min_work_it);
                // Now make the chunk if it's ours
                if (assigned_process == ProcessorRank)
                {
                    chunks.emplace_back(config_index, config_index + current_num_configs, csf_start,
                                        current_num_rows, Nsmall);
                }

                processor_work_sizes[assigned_process] += current_chunk_work_units;
            }

            config_index += current_num_configs;
            csf_start += current_num_rows;
            most_chunk_rows = mmax(most_chunk_rows, current_num_rows);
        }

        // Print some diagnostics about the workload balancing
        // What is this process's chunk workload?
        *logstream << "\nThis process has " << processor_work_sizes[ProcessorRank] << " work units and "
                   << num_big_chunks
                   << " big chunks.\n";

        // How big is the workload imbalance?
        auto min_work_it = std::min_element(processor_work_sizes.begin(), processor_work_sizes.end());
        auto max_work_it = std::max_element(processor_work_sizes.begin(), processor_work_sizes.end());
        double imbalance = 100.0 * ((double) (*max_work_it) - (double) (*min_work_it)) / ((double) (*max_work_it));

        *logstream << "Minimum workload: " << *min_work_it;
        *logstream << "\nMaximum workload: " << *max_work_it;
        *logstream << "\nThe relative workload imbalance across MPI processes is " << imbalance << "%" << std::endl;
    }

    // Loop through my chunks
    auto config_it = configs->begin();
    RelativisticConfigList::const_iterator configsubsetend_it = configs->small_end();
    unsigned int configsubsetend = configs->small_size();

#ifdef AMBIT_USE_OPENMP
    Eigen::setNbThreads(1);
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
            auto config_j_index = 0;
            RelativisticConfigList::const_iterator config_jend;
            if(config_index < configsubsetend)
                config_jend = config_it;
            else
                config_jend = configsubsetend_it;

            RowMajorMatrix H_proj;

#ifdef AMBIT_USE_OPENMP
            #pragma omp task untied \
                             default(none) \
                             shared(M) \
                             firstprivate(leading_config_i, config_it, config_jt, config_jend, \
                                          current_chunk, config_index, config_j_index, H_proj)
            {
#endif
            // All configs up to the diagonal (or smallsize)
            while(config_jt != config_jend)
            {
                bool leading_config_j = H_three_body && std::binary_search(leading_configs->first.begin(), leading_configs->first.end(), NonRelConfiguration(*config_jt));

                int config_diff_num = config_it->GetConfigDifferencesCount(*config_jt);
                bool do_three_body = (leading_config_i || leading_config_j) && (config_diff_num <= 3);

                // Check that the number of differences is small enough
                if(do_three_body || (config_diff_num <= 2))
                {
                    auto start_time = std::chrono::high_resolution_clock::now();

                    H_proj = RowMajorMatrix::Zero(config_it->projection_size(), config_jt->projection_size());

                    // Loop through projections
                    auto proj_it = config_it.projection_begin();
                    int pi = 0;

                    while(proj_it != config_it.projection_end())
                    {
                        auto proj_jt = config_jt.projection_begin();
                        int pj = 0;

                        while(proj_jt != config_jt.projection_end())
                        {
                            if(do_three_body)
                                H_proj(pi, pj) = H_three_body->GetMatrixElement(*proj_it, *proj_jt);
                            else
                                H_proj(pi, pj) = H_two_body->GetMatrixElement(*proj_it, *proj_jt);

                            proj_jt++; pj++;
                        }
                        proj_it++; pi++;
                    }

                    // Get CSF data
                    Eigen::Map<const RowMajorMatrix> angular_i_mapped(config_it->GetCSFs(), config_it->projection_size(), config_it->NumCSFs());
                    Eigen::Map<const RowMajorMatrix> angular_j_mapped(config_jt->GetCSFs(), config_jt->projection_size(), config_jt->NumCSFs());

                    M.block(config_it.csf_offset()-current_chunk.start_row, config_jt.csf_offset(), config_it->NumCSFs(), config_jt->NumCSFs())
                        = angular_i_mapped.transpose() * H_proj * angular_j_mapped;
                }
                config_jt++;
                config_j_index++;
            } // while (config_jt)
#ifdef AMBIT_USE_OPENMP
            } // OMP task
#endif

#ifdef AMBIT_USE_OPENMP
            #pragma omp task untied \
                             default(none) \
                             shared(M, D) \
                             firstprivate(current_chunk, config_it, config_index, leading_config_i, configsubsetend, H_proj)
            {
#endif
            auto start_time = std::chrono::high_resolution_clock::now();

            // Create matrix <proj_i | H | proj_j>
            H_proj = RowMajorMatrix::Zero(config_it->projection_size(), config_it->projection_size());

            // Loop through projections to get upper part of matrix H_proj
            auto proj_it = config_it.projection_begin();
            int pi = 0;
            while(proj_it != config_it.projection_end())
            {
                auto proj_jt = proj_it;
                int pj = pi;
                while (proj_jt != config_it.projection_end())
                {
                    if (leading_config_i)
                        H_proj(pi, pj) = H_three_body->GetMatrixElement(*proj_it, *proj_jt);
                    else
                        H_proj(pi, pj) = H_two_body->GetMatrixElement(*proj_it, *proj_jt);

                    proj_jt++; pj++;
                }
                proj_it++; pi++;
            }

            // Get CSF data
            Eigen::Map<const RowMajorMatrix> angular_mapped(config_it->GetCSFs(), config_it->projection_size(), config_it->NumCSFs());

            if(config_index < configsubsetend)
            {
                M.block(config_it.csf_offset()-current_chunk.start_row, config_it.csf_offset(), config_it->NumCSFs(), config_it->NumCSFs())
                    = angular_mapped.transpose() * H_proj.selfadjointView<Eigen::Upper>() * angular_mapped;
            }
            else
            {
                int diag_offset = current_chunk.start_row + current_chunk.num_rows - current_chunk.diagonal.rows();
                D.block(config_it.csf_offset()-diag_offset, config_it.csf_offset()-diag_offset, config_it->NumCSFs(), config_it->NumCSFs())
                    = angular_mapped.transpose() * H_proj.selfadjointView<Eigen::Upper>() * angular_mapped;
            }
#ifdef AMBIT_USE_OPENMP
            } // OMP diagonal task
#endif
            config_it++;
        } // Configs in chunk
    } // Chunks loop
#ifdef AMBIT_USE_OPENMP
    } // OMP single
    }  // OMP Parallel
    Eigen::setNbThreads(0);
#endif

    std::sort(chunks.begin(), chunks.end(), [](const MatrixChunk& left, const MatrixChunk& right) -> bool
        {   return (left.start_row < right.start_row); });

    for(auto& matrix_section: chunks)
        matrix_section.Symmetrize();
}

LevelVector HamiltonianMatrix::SolveMatrix(pHamiltonianID hID, unsigned int num_solutions, const std::string& filename)
{
    LevelVector levelvec(hID, configs);
    unsigned int NumSolutions = mmin(num_solutions, N);

    if(NumSolutions == 0)
    {
        *outstream << "\nNo solutions" << std::endl;
    }
    else if(N <= SMALL_MATRIX_LIM || NumSolutions > MANY_LEVELS_LIM)
    {
        // Use Eigen or SCALAPACK
        RowMajorMatrix M;
        RowMajorMatrix* pM;

        if(N <= SMALL_MATRIX_LIM)
        {   // There is only a single chunk: no need to copy
            *outstream << "; Finding solutions using Eigen..." << std::endl;
            pM = &chunks.front().chunk;
        }
        else if(NumProcessors == 1)
        {   // Copy all chunks to a single matrix.
           *outstream << "; Attempting to reallocate matrix and find solutions using Eigen..." << std::endl;

            M = RowMajorMatrix::Zero(N, N);
            pM = &M;

            for(auto& chunk: chunks)
            {
                M.block(chunk.start_row, 0, chunk.chunk.rows(), chunk.chunk.cols()) = chunk.chunk;
                if(chunk.diagonal.size())
                {
                    int diag_offset = chunk.start_row + chunk.num_rows - chunk.diagonal.rows();
                    M.block(diag_offset, diag_offset, chunk.diagonal.rows(), chunk.diagonal.cols()) = chunk.diagonal;
                }
            }
        }
        else
        {   // Write matrix file, clear current Hamiltonian to make space,
            if(filename.empty())
                return levelvec;
            Write(filename);
            Clear();

#ifdef AMBIT_USE_SCALAPACK
            // Read and solve using ScalapackMatrix
            return SolveMatrixScalapack(hID, num_solutions, filename);
#else
            // Read to single processor
            *outstream << "; Attempting to reallocate matrix and find solutions using Eigen..." << std::endl;
            Read(filename, configs->size(), true);
            pM = &chunks.front().chunk;
#endif
        }
        levelvec.Resize(NumSolutions, N);

        if(ProcessorRank == 0)
        {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(*pM);
            const Eigen::VectorXd& E = es.eigenvalues();
            const Eigen::MatrixXd& V = es.eigenvectors();

            for(unsigned int i = 0; i < NumSolutions; i++)
            {
                levelvec.eigenvalues[i] = E(i);
                levelvec.eigenvectors.row(i) = V.col(i);
            }
        }
#ifdef AMBIT_USE_MPI
        if(NumProcessors > 1)
        {   // broadcast results
            MPI_Bcast(levelvec.eigenvalues.data(), NumSolutions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(levelvec.eigenvectors.data(), NumSolutions * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
#endif
    }
    else
    {   *outstream << "; Finding solutions using Davidson..." << std::endl;
        levelvec.Resize(NumSolutions, N);

        double* V = levelvec.eigenvectors.data();
        double* E = levelvec.eigenvalues.data();

        Eigensolver solver;
        #ifdef AMBIT_USE_MPI
            solver.MPISolveLargeSymmetric(this, E, V, N, NumSolutions);
        #else
            solver.SolveLargeSymmetric(this, E, V, N, NumSolutions);
        #endif
    }

    return levelvec;
}

#ifdef AMBIT_USE_SCALAPACK
LevelVector HamiltonianMatrix::SolveMatrixScalapack(pHamiltonianID hID, unsigned int num_solutions, const std::string& filename, std::optional<double> energy_limit)
{
    LevelVector levelvec(hID, configs);
    unsigned int NumSolutions = mmin(num_solutions, N);

    if(NumSolutions == 0)
    {
        *outstream << "\nNo solutions" << std::endl;
    }
    else
    {   *outstream << "; Finding solutions using ScaLAPACK ..." << std::endl;
        Clear();

        ScalapackMatrix SM(N);
        *logstream << "ScalapackMatrix allocated..." << std::flush;
        SM.ReadLowerTriangle(filename);
        *logstream << " read..." << std::flush;

        // Diagonalise
        levelvec.eigenvalues.resize(N); // All eigenvalues
        double* E = levelvec.eigenvalues.data();
        SM.Diagonalise(E);
        *logstream << " diagonalised." << std::endl;

        // Cut off num_solutions
        if(energy_limit)
            for(int i = 0; i < NumSolutions; i++)
            {
                if(E[i] > energy_limit)
                {
                    NumSolutions = i;
                    break;
                }
            };

        // Get levels
        levelvec.Resize(NumSolutions, N);

        // Do a few sets of eigenvectors at a time, order 10M elements
        // Scalapack returns column-major vectors
        unsigned int column_begin = 0;
        unsigned int num_columns_per_step = mmax(10000000/N, 1);

        ColMajorMatrix V(N, num_columns_per_step);

        while(column_begin < NumSolutions)
        {
            unsigned int column_end = mmin(column_begin + num_columns_per_step, NumSolutions);
            SM.GetColumns(column_begin, column_end, V.data());
            unsigned int num_eigenvectors = column_end-column_begin;

            levelvec.eigenvectors.middleRows(column_begin, num_eigenvectors) = V.leftCols(num_eigenvectors).transpose();
            column_begin = column_end;
        }
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

bool HamiltonianMatrix::Read(const std::string& filename, unsigned int configs_per_chunk, bool only_root_process)
{
    FILE* fp = file_err_handler->fopen(filename.c_str(), "rb");
    if(!fp)
        return false;

    // NOTE: The processes don't yet know which chunks they "own", nor which rows that corresponds
    // to because that's all calculated in HamiltonianMatrix::GenerateMatrix() and this function is
    // intended to help us avoid having to call that. We'll need to re-calculate the chunk
    // distribution which may or may not be different from what we used in the run that generated
    // the `.matrix` file (e.g. if we wrote the file using N processes and are reading with M != N
    // processes). Fortunately, if we're reading a complete matrix, then we clearly don't need to
    // worry about workload imbalances in GenerateMatrix() (it's already been generated!), so we
    // can just do a really simply round-robin assignment of chunks to processes.
    if(chunks.empty())
    {
        if((N <= SMALL_MATRIX_LIM) || only_root_process)
        {
            configs_per_chunk = configs->size();
        }

        // Total number of chunks = ceiling(number of configs/configs_per_chunk)
        unsigned int total_num_chunks = (configs->size() + configs_per_chunk - 1)/configs_per_chunk;
        // Now loop through and construct this rank's chunks (i.e. allocate the space but don't
        // actually read anything yet)
        auto config_it = configs->begin();
        unsigned int config_index = 0;
        unsigned int csf_start = 0;

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
            // Now make the chunk if it's ours
            if(chunk_index%NumProcessors == ProcessorRank)
            {
                chunks.emplace_back(config_index, config_index+current_num_configs, csf_start,
                                    current_num_rows, Nsmall);
            }
            config_index += current_num_configs;
            csf_start += current_num_rows;
            most_chunk_rows = mmax(most_chunk_rows, current_num_rows);
        }
    }

    // Read size of matrix
    unsigned int matrix_N;
    file_err_handler->fread(&matrix_N, sizeof(unsigned int), 1, fp);

    if(matrix_N != N)
    {   *errstream << "HamiltonianMatrix::Read: matrix has incorrect dimension (number of CSFs)." << std::endl;
        *errstream << "Expected " << N << " but got " << matrix_N << std::endl;
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

        long long int row = 0;
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
                int err;
                err = MPI_Recv(&buf, Nsmall*most_chunk_rows, MPI_DOUBLE, MPI_ANY_SOURCE, row, MPI_COMM_WORLD, &status);
                if(err)
                    *errstream << "HamiltonianMatrix::Write: MPI_Error " << err << std::endl;

                // Get number of rows in chunk
                int data_count;
                MPI_Get_count(&status, MPI_DOUBLE, &data_count);
                if(data_count/Nsmall > int(Nsmall) - row)
                    num_rows = data_count/Nsmall;
                else
                    num_rows = (-row + sqrtl(row * row + 4LL * (long long int)(data_count)))/2;

                if(num_rows * mmin(row + num_rows, Nsmall) != data_count)
                {   *errstream << "HamiltonianMatrix::Write: received incorrect chunk size." << std::endl;
                    *errstream << "  row=" << row << " num_rows=" << num_rows << " Nsmall=" << Nsmall << " data_count=" << data_count << std::endl;
                    exit(1);
                }

                pbuf = buf;

                // Receive diagonal
                if(row + num_rows > Nsmall)
                {
                    MPI_Recv(&diagbuf, most_chunk_rows*most_chunk_rows, MPI_DOUBLE, MPI_ANY_SOURCE, row+N, MPI_COMM_WORLD, &status);

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
                    MPI_Send(chunk_it->diagonal.data(), chunk_it->diagonal.size(), MPI_DOUBLE, 0, row+N, MPI_COMM_WORLD);

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
