#include "Include.h"
#include "HamiltonianMatrix.h"
#include "Universal/SmallMatrix.h"
#include "Universal/SymMatrix.h"
#include "HartreeFock/SingleParticleWavefunction.h"
#include "Universal/Eigensolver.h"
#include "Universal/MathConstant.h"
#include "GFactor.h"

// Don't bother with davidson method if smaller than this limit
#define SMALL_MATRIX_LIM 1

// Include this define for the box diagrams of "wrong" parity.
//#define INCLUDE_EXTRA_BOX_DIAGRAMS

// Include this define to only include sigma3 when both states are leading configurations
// (instead of just one).
//#define SIGMA3_AND

HamiltonianMatrix::HamiltonianMatrix(pOneElectronIntegrals hf, pTwoElectronCoulombOperator coulomb, pRelativisticConfigListConst relconfigs):
    H_two_body(nullptr), configs(relconfigs) //, include_sigma3(false)
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
    }

    // Loop through my chunks
    for(auto& current_chunk: chunks)
    {
        Eigen::MatrixXd& M = current_chunk.chunk;

        // Loop through configs for this chunk
#ifdef AMBIT_USE_OPENMP
        #pragma omp parallel for private(config_it)
#endif
        for(unsigned int config_index = current_chunk.config_indices.first; config_index < current_chunk.config_indices.second; config_index++)
        {
            auto current_config_pair = (*configs)[config_index];
            config_it = current_config_pair.first;
            int csf_offset_i = current_config_pair.second;

            // Loop through projections
            auto proj_it = config_it->projection_begin();
            while(proj_it != config_it->projection_end())
            {
                // Loop through the rest of the configs
                auto config_jt = config_it;
                int csf_offset_j = csf_offset_i;
                while(config_jt != configs->end())
                {
                    RelativisticConfiguration::const_projection_iterator proj_jt;
                    if(config_jt == config_it)
                        proj_jt = proj_it;
                    else
                        proj_jt = config_jt->projection_begin();

                    while(proj_jt != config_jt->projection_end())
                    {
                        double operatorH = H_two_body->GetMatrixElement(*proj_it, *proj_jt);

                        if(fabs(operatorH) > 1.e-15)
                        {
                            for(auto coeff_i = proj_it.CSF_begin(csf_offset_i); coeff_i != proj_it.CSF_end(csf_offset_i); coeff_i++)
                            {
                                RelativisticConfigList::const_CSF_iterator start_j = proj_jt.CSF_begin(csf_offset_j);

                                if(proj_it == proj_jt)
                                    start_j = coeff_i;

                                for(auto coeff_j = start_j; coeff_j != proj_jt.CSF_end(csf_offset_j); coeff_j++)
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

                    csf_offset_j += config_jt->NumCSFs();
                    config_jt++;
                }

                proj_it++;
            }
        } // Configs in chunk
    } // Chunks

    for(auto& matrix_section: chunks)
        matrix_section.Symmetrize();
}

void HamiltonianMatrix::SolveMatrix(const Symmetry& sym, unsigned int num_solutions, pLevelMap levels)
{
    if(N == 0)
    {   *outstream << "\nNo solutions" << std::endl;
        return;
    }
    
    unsigned int NumSolutions = mmin(num_solutions, N);
    
    *outstream << "; Finding solutions..." << std::endl;
    
    if(N <= SMALL_MATRIX_LIM)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(chunks.front().chunk);
        const Eigen::VectorXd& E = es.eigenvalues();
        const Eigen::MatrixXd& V = es.eigenvectors();
        
        for(unsigned int i = 0; i < NumSolutions; i++)
        {
            pLevel level(new Level(E(i), V.col(i).data(), configs, N));
            (*levels)[LevelID(sym, i)] = level;
        }
    }
    else
    {   // Using Davidson method
        double* V = new double[NumSolutions * N];
        double* E = new double[NumSolutions];
        
        Eigensolver solver;
        #ifdef _MPI
            solver.MPISolveLargeSymmetric(this, E, V, N, NumSolutions);
        #else
            solver.SolveLargeSymmetric(this, E, V, N, NumSolutions);
        #endif

        for(unsigned int i = 0; i < NumSolutions; i++)
        {
            pLevel level(new Level(E[i], (V + N * i), configs, N));
            (*levels)[LevelID(sym, i)] = level;
        }
        
        delete[] E;
        delete[] V;
    }
}

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
    return;
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
