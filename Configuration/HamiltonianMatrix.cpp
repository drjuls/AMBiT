#include "Include.h"
#include "HamiltonianMatrix.h"
#include "HartreeFock/Orbital.h"
#include "Universal/Eigensolver.h"
#include "Universal/MathConstant.h"
#include "Universal/ScalapackMatrix.h"
#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif

// Don't bother with davidson method if smaller than this limit
#define SMALL_MATRIX_LIM 200

// Don't bother with davidson method if number of solutions requested is larger than this
// and ScaLAPACK is available.
#define MANY_LEVELS_LIM   50

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
        *logstream << " " << N << "x" << Nsmall << std::flush;
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
    unsigned int num_chunks = (configs->size() + configs_per_chunk - 1)/configs_per_chunk;

    // Divide up chunks
    auto config_it = configs->begin();
    unsigned int config_index = 0;
    unsigned int csf_start = 0;
    for(int chunk_index = 0; chunk_index < num_chunks; chunk_index++)
    {
        // Get chunk num_rows and number of configs
        unsigned int current_num_rows = 0;
        unsigned int current_num_configs = 0;
        while(config_it != configs->end() && current_num_configs < configs_per_chunk)
        {
            current_num_rows += config_it->NumCSFs();
            current_num_configs++;
            config_it++;
        }

        if(current_num_rows == 0)
            break;

        // Make chunk
        if(chunk_index%NumProcessors == ProcessorRank)
            chunks.emplace_back(config_index, config_index+current_num_configs, csf_start, current_num_rows, Nsmall);

        config_index += current_num_configs;
        csf_start += current_num_rows;
        most_chunk_rows = mmax(most_chunk_rows, current_num_rows);
    }

    // Loop through my chunks
    RelativisticConfigList::const_iterator configsubsetend_it = configs->small_end();
    unsigned int configsubsetend = configs->small_size();

    for(auto& current_chunk: chunks)
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& M = current_chunk.chunk;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& D = current_chunk.diagonal;

        // Loop through configs for this chunk
#ifdef AMBIT_USE_OPENMP
        #pragma omp parallel for private(config_it)
#endif
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
                                operatorH = H_three_body->GetMatrixElement(*proj_it, *proj_jt);
                            else
                                operatorH = H_two_body->GetMatrixElement(*proj_it, *proj_jt);

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
            }

            // Diagonal
            if(config_index >= configs->small_size())
            {
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
                }
            }

            config_it++;
        } // Configs in chunk
    } // Chunks

    for(auto& matrix_section: chunks)
        matrix_section.Symmetrize();
}

LevelVector HamiltonianMatrix::SolveMatrix(pHamiltonianID hID, unsigned int num_solutions)
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
    #ifdef AMBIT_USE_SCALAPACK
        else if(NumSolutions > MANY_LEVELS_LIM && Nsmall == N)
        {
            levelvec = SolveMatrixScalapack(hID, NumSolutions, false);
        }
    #endif
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

void HamiltonianMatrix::Write(const std::string& filename) const
{
    FILE* fp;

    auto chunk_it = chunks.begin();

    // Send rows to root node, which writes them sequentially.
    if(ProcessorRank == 0)
    {
        // Write size of matrix.
        fp = fopen(filename.c_str(), "wb");
        fwrite(&N, sizeof(unsigned int), 1, fp);
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
                fwrite(pbuf, sizeof(double), mmin(row + i + 1, Nsmall), fp);
                pbuf += mmin(row + num_rows, Nsmall);   // Move to next row

                // Write diagonal
                if(diag_rows && row + i >= Nsmall)
                {
                    if(gap)
                        fwrite(zeros.data(), sizeof(double), gap, fp);
                    fwrite(pdiag, sizeof(double), row + i + 1 - Nsmall - gap, fp);
                    pdiag += diag_rows;
                }
            }

            row += num_rows;
        }

    #ifdef AMBIT_USE_MPI
        // Send finished signal
        MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);
    #endif

        fclose(fp);
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
