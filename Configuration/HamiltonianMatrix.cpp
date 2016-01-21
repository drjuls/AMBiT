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

// Include this define for the box diagrams of "wrong" parity.
//#define INCLUDE_EXTRA_BOX_DIAGRAMS

// Include this define to only include sigma3 when both states are leading configurations
// (instead of just one).
//#define SIGMA3_AND

HamiltonianMatrix::HamiltonianMatrix(pHFIntegrals hf, pTwoElectronCoulombOperator coulomb, pRelativisticConfigList relconfigs):
    H_two_body(nullptr), configs(relconfigs), most_chunk_rows(0) //, include_sigma3(false)
{
    // Set up Hamiltonian operator
    H_two_body = pTwoBodyHamiltonianOperator(new TwoBodyHamiltonianOperator(hf, coulomb));

    // Set up matrix
    N = configs->NumCSFs();

    *logstream << " " << N << " " << std::flush;
    *outstream << " Number of CSFs = " << N << std::flush;
}

HamiltonianMatrix::~HamiltonianMatrix()
{}

void HamiltonianMatrix::GenerateMatrix()
{
    chunks.clear();

    unsigned int configs_per_chunk = 4;

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
            chunks.emplace_back(config_index, config_index+current_num_configs, csf_start, current_num_rows, N);

        config_index += current_num_configs;
        csf_start += current_num_rows;
        most_chunk_rows = mmax(most_chunk_rows, current_num_rows);
    }

    // Loop through my chunks
    for(auto& current_chunk: chunks)
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& M = current_chunk.chunk;

        // Loop through configs for this chunk
#ifdef AMBIT_USE_OPENMP
        #pragma omp parallel for private(config_it)
#endif
        config_it = (*configs)[current_chunk.config_indices.first];
        for(unsigned int config_index = current_chunk.config_indices.first; config_index < current_chunk.config_indices.second; config_index++)
        {
            // Loop through the rest of the configs
            auto config_jt = config_it;
            while(config_jt != configs->end())
            {
                // Check that the number of differences is small enough
                if(config_it->GetConfigDifferencesCount(*config_jt) <= 2)
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
                                        int i = coeff_i.index() - current_chunk.start_row;
                                        int j = coeff_j.index() - current_chunk.start_row;

                                        if(i < j)
                                            M(i, j) += operatorH * (*coeff_i) * (*coeff_j);
                                        else if(i > j)
                                            M(j, i) += operatorH * (*coeff_i) * (*coeff_j);
                                        else if(proj_it == proj_jt)
                                            M(i, j) += operatorH * (*coeff_i) * (*coeff_j);
                                        else
                                            M(i, j) += 2. * operatorH * (*coeff_i) * (*coeff_j);
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
            config_it++;
        } // Configs in chunk
    } // Chunks

    for(auto& matrix_section: chunks)
        matrix_section.Symmetrize();
}

LevelVector HamiltonianMatrix::SolveMatrix(pHamiltonianID hID, unsigned int num_solutions)
{
    LevelVector levels;
    if(hID->GetRelativisticConfigList() == nullptr)
        hID->SetRelativisticConfigList(configs);

    unsigned int NumSolutions = mmin(num_solutions, N);

    if(NumSolutions == 0)
    {
        *outstream << "\nNo solutions" << std::endl;
    }
    else
    {   if(N <= SMALL_MATRIX_LIM && NumProcessors == 1)
        {
            *outstream << "; Finding solutions using Eigen..." << std::endl;
            levels.reserve(NumSolutions);

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(chunks.front().chunk);
            const Eigen::VectorXd& E = es.eigenvalues();
            const Eigen::MatrixXd& V = es.eigenvectors();

            for(unsigned int i = 0; i < NumSolutions; i++)
            {
                levels.push_back(std::make_shared<Level>(E(i), V.col(i).data(), hID, N));
            }
        }
    #ifdef AMBIT_USE_SCALAPACK
        else if(NumSolutions > MANY_LEVELS_LIM)
        {
            levels = SolveMatrixScalapack(hID, NumSolutions, false);
        }
    #endif
        else
        {   *outstream << "; Finding solutions using Davidson..." << std::endl;
            levels.reserve(NumSolutions);

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
                levels.push_back(std::make_shared<Level>(E[i], (V + N * i), hID, N));
            }

            delete[] E;
            delete[] V;
        }
    }

    return levels;
}

#ifdef AMBIT_USE_SCALAPACK
LevelVector HamiltonianMatrix::SolveMatrixScalapack(pHamiltonianID hID, unsigned int num_solutions, bool use_energy_limit, double energy_limit)
{
    LevelVector levels;
    if(hID->GetRelativisticConfigList() == nullptr)
        hID->SetRelativisticConfigList(configs);

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
        SM.ReadTriangle(filename);

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
        levels.reserve(NumSolutions);
        unsigned int column_begin = 0;
        unsigned int num_columns_per_step = 1000000/N;
        double* V = new double[N * num_columns_per_step];  // Eigenvectors

        while(column_begin < NumSolutions)
        {
            unsigned int column_end = mmin(column_begin + num_columns_per_step, NumSolutions);
            SM.GetColumns(column_begin, column_end, V);

            double* pV = V;
            while(column_begin < column_end)
            {   levels.push_back(std::make_shared<Level>(E[column_begin], pV, hID, N));
                column_begin++;
                pV += N;
            }
        }

        delete[] E;
        delete[] V;
    }

    return levels;
}
#endif

std::ostream& operator<<(std::ostream& stream, const HamiltonianMatrix& matrix)
{
    for(auto& matrix_section: matrix.chunks)
    {
        // Each row separately
        for(unsigned int row = 0; row < matrix_section.num_rows; row++)
        {
            // Leading zeros
            stream << Eigen::VectorXd::Zero(matrix_section.start_row + row).transpose() << " ";

            // Upper triangular matrix part of row
            stream << matrix_section.chunk.block(row, row, 1, matrix_section.chunk.cols() - row) << "\n";
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

    #ifdef AMBIT_USE_MPI
        double buf[N * most_chunk_rows];
    #endif

        unsigned int row = 0;
        while(row < N)
        {
            int num_rows = 0;

            // My chunk!
            if(row == chunk_it->start_row)
            {
                pbuf = chunk_it->chunk.data();
                num_rows = chunk_it->num_rows;
                chunk_it++;
            }
        #ifdef AMBIT_USE_MPI
            else
            {   // Broadcast row number
                MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);

                // Receive chunk
                MPI_Status status;
                MPI_Recv(&buf, N*most_chunk_rows, MPI_DOUBLE, MPI_ANY_SOURCE, row, MPI_COMM_WORLD, &status);

                // Get number of rows in chunk
                MPI_Get_count(&status, MPI_DOUBLE, &num_rows);
                if(num_rows%(N-row))
                    *errstream << "HamiltonianMatrix::Write: received chunk size not a multiple of (N - row)." << std::endl;
                num_rows = num_rows/(N-row);

                pbuf = buf;
            }
        #endif

            for(int i = 0; i < num_rows; i++)
            {
                fwrite(pbuf, sizeof(double), N - row - i, fp);
                pbuf += (N - row + 1);  // Move one more column and one more row along
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
                MPI_Send(chunk_it->chunk.data(), (N - row) * chunk_it->num_rows, MPI_DOUBLE, 0, row, MPI_COMM_WORLD);
                chunk_it++;
            }
        }
    #endif
    }
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

        // Upper triangular part
        c_mapped.middleRows(start, matrix_section.num_rows).noalias()
            += matrix_section.chunk * b_mapped.bottomRows(N-start);

        // Lower triangular part
        unsigned int lower_start = start + matrix_section.num_rows;
        if(lower_start < N)
        {   c_mapped.bottomRows(N - lower_start).noalias()
                += matrix_section.chunk.rightCols(N - lower_start).transpose() * b_mapped.middleRows(start, matrix_section.num_rows);
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
        diag_mapped.segment(matrix_section.start_row, matrix_section.num_rows).noalias()
            = matrix_section.chunk.diagonal();
    }
}

//double HamiltonianMatrix::GetSigma3(const Projection& first, const Projection& second) const
//{
//    const std::set<Configuration>* leading_configs = confgen->GetLeadingConfigs();
//
//#ifdef SIGMA3_AND
//    // Check that first AND second are leading configurations
//    if((leading_configs->find(first.GetNonRelConfiguration()) == leading_configs->end()) ||
//       (leading_configs->find(second.GetNonRelConfiguration()) == leading_configs->end()))
//        return 0.;
//#else
//    // Check that first OR second is a leading configuration
//    if((leading_configs->find(first.GetNonRelConfiguration()) == leading_configs->end()) &&
//       (leading_configs->find(second.GetNonRelConfiguration()) == leading_configs->end()))
//        return 0.;
//#endif
//
//    unsigned int diff[6];
//    int numdiff = Projection::GetProjectionDifferences3(first, second, diff);
//
//    int sign;
//    if(numdiff >= 0)
//        sign = 1;
//    else
//        sign = -1;
//
//    double value = 0.;
//
//    if(numdiff == 0)
//    {
//        // Sum(i < j < k) Sigma3(ijk, ijk)
//        for(unsigned int i=0; i<first.size(); i++)
//        {
//            for(unsigned int j=i+1; j<first.size(); j++)
//            {
//                for(unsigned int k=j+1; k<first.size(); k++)
//                {
//                    value += Sigma3(first[i], first[j], first[k], first[i], first[j], first[k]);
//                }
//            }
//        }
//    }
//    else if(abs(numdiff) == 1)
//    {
//        const ElectronInfo& f1 = first[diff[0]];
//        const ElectronInfo& s1 = second[diff[1]];
//
//        // Sum(i < j) Sigma3(aij, bij)
//        for(unsigned int i=0; i<first.size(); i++)
//        {
//            for(unsigned int j=i+1; j<first.size(); j++)
//            {
//                if((i != diff[0]) && (j != diff[0]))
//                {
//                    value += sign * Sigma3(f1, first[i], first[j], s1, first[i], first[j]);
//                }
//            }
//        }
//    }
//    else if(abs(numdiff) == 2)
//    {
//        const ElectronInfo& f1 = first[diff[0]];
//        const ElectronInfo& s1 = second[diff[1]];
//        const ElectronInfo& f2 = first[diff[2]];
//        const ElectronInfo& s2 = second[diff[3]];
//
//        // Sum(i) Sigma3(abi, cdi)
//        for(unsigned int i=0; i<first.size(); i++)
//        {
//            if((i != diff[0]) && (i != diff[2]))
//               value += sign * Sigma3(f1, f2, first[i], s1, s2, first[i]);
//        }
//    }
//    else if(abs(numdiff) == 3)
//    {
//        const ElectronInfo& f1 = first[diff[0]];
//        const ElectronInfo& s1 = second[diff[1]];
//        const ElectronInfo& f2 = first[diff[2]];
//        const ElectronInfo& s2 = second[diff[3]];
//        const ElectronInfo& f3 = first[diff[4]];
//        const ElectronInfo& s3 = second[diff[5]];
//
//        // Sigma3(abc, def)
//        value = sign * Sigma3(f1, f2, f3, s1, s2, s3);
//    }
//
//    return value;
//}
//
//double HamiltonianMatrix::Sigma3(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
//           const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const
//{
//    // Check momentum projections
//    if(e1.TwoM() + e2.TwoM() + e3.TwoM() != e4.TwoM() + e5.TwoM() + e6.TwoM())
//        return 0.;
//
//    // Check parity
//    if((e1.L() + e2.L() + e3.L() + e4.L() + e5.L() + e6.L())%2)
//        return 0.;
//
//    double value = 0.;
//
//    // The sign changes for odd permutations
//    value =   Sigma3LinePermutations(e1, e2, e3, e4, e5, e6)
//            + Sigma3LinePermutations(e1, e2, e3, e5, e6, e4)
//            + Sigma3LinePermutations(e1, e2, e3, e6, e4, e5)
//            - Sigma3LinePermutations(e1, e2, e3, e5, e4, e6)
//            - Sigma3LinePermutations(e1, e2, e3, e6, e5, e4)
//            - Sigma3LinePermutations(e1, e2, e3, e4, e6, e5);
//
//    return value;
//}
//
///** This function does the line permutations, putting the pairs on different levels
//    of the three-body interaction.
//    */
//inline double HamiltonianMatrix::Sigma3LinePermutations(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
//           const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const
//{
//    double value = 0.;
//
//    // There are no sign changes, since there are the same number
//    // of permutations on both sides
//
//    value =   sigma3calc->GetSecondOrderSigma3(e1, e2, e3, e4, e5, e6)
//            + sigma3calc->GetSecondOrderSigma3(e2, e3, e1, e5, e6, e4)
//            + sigma3calc->GetSecondOrderSigma3(e3, e1, e2, e6, e4, e5)
//            + sigma3calc->GetSecondOrderSigma3(e3, e2, e1, e6, e5, e4)
//            + sigma3calc->GetSecondOrderSigma3(e2, e1, e3, e5, e4, e6)
//            + sigma3calc->GetSecondOrderSigma3(e1, e3, e2, e4, e6, e5);
//
//    return value;
//}
