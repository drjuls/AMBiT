#include "Include.h"
#include "HamiltonianMatrix.h"
#include "Universal/SmallMatrix.h"
#include "Universal/SymMatrix.h"
#include "HartreeFock/SingleParticleWavefunction.h"
#include "Universal/Eigensolver.h"
#include "Universal/MathConstant.h"
//#include "ConfigFileGenerator.h"
#include "GFactor.h"

#define SMALL_MATRIX_LIM 2000

// Include this define for the box diagrams of "wrong" parity.
//#define INCLUDE_EXTRA_BOX_DIAGRAMS

// Include this define to only include sigma3 when both states are leading configurations
// (instead of just one).
//#define SIGMA3_AND

HamiltonianMatrix::HamiltonianMatrix(pHFElectronOperator hf, const CIIntegrals& coulomb_integrals, pRelativisticConfigListConst relconfigs):
    H_two_body(nullptr), configs(relconfigs), M(nullptr) //, include_sigma3(false)
{
    // Set up Hamiltonian operator
    pTwoElectronCoulombOperator coulomb(new TwoElectronCoulombOperator(coulomb_integrals));
    H_two_body = pTwoBodyHamiltonianOperator(new TwoBodyHamiltonianOperator(hf, coulomb));

    // Set up matrix
    N = configs->NumCSFs();

    *logstream << " " << N << " " << std::flush;
    *outstream << " Number of CSFs = " << N << std::endl;
}

HamiltonianMatrix::~HamiltonianMatrix()
{
    if(M)
        delete M;
}

void HamiltonianMatrix::GenerateMatrix()
{
    if(M == nullptr)
    {   if(N <= SMALL_MATRIX_LIM)
            M = new SmallMatrix(N);
        else
            M = new SymMatrix(N);
    }
    else
        M->Clear();

    M->WriteMode(true);

    // Loop through projections
    auto proj_it = configs->projection_begin();
    while(proj_it != configs->projection_end())
    {
        auto proj_jt = proj_it;

        while(proj_jt != configs->projection_end())
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

                        if(i < j)
                            M->At(i, j) += operatorH * (*coeff_i) * (*coeff_j);
                        else if(i > j)
                            M->At(j, i) += operatorH * (*coeff_i) * (*coeff_j);
                        else if(proj_it == proj_jt)
                            M->At(i, j) += operatorH * (*coeff_i) * (*coeff_j);
                        else
                            M->At(i, j) += 2. * operatorH * (*coeff_i) * (*coeff_j);
                    }
                }
            }
            proj_jt++;
        }
        proj_it++;
    }
}

void HamiltonianMatrix::WriteToFile(const std::string& filename)
{
    return;
}

void HamiltonianMatrix::PollMatrix()
{
    M->WriteMode(false);

    unsigned int i, j, count;
    double value;
    unsigned int range[10];
    for(i=0; i<10; i++)
        range[i] = 0;

    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
        {
            value = fabs(M->At(i, j));
            if(value >= 1.)
                range[9]++;
            else
            {   count = 0;
                while(value > 1.e-16)
                {   value = value/100.;
                    count++;
                }
                range[count]++;
            }
        }

    for(i=0; i<10; i++)
        *outstream << i << " " << range[i] << " " << double(range[i])/double(N*N)*100. << std::endl;
}

void HamiltonianMatrix::SolveMatrix(const Symmetry& sym, unsigned int num_solutions, pLevelMap levels)
{
    if(N == 0)
    {   *outstream << "\nNo solutions" << std::endl;
        return;
    }

    M->WriteMode(false);

    *outstream << "\nFinding solutions" << std::endl;

    unsigned int NumSolutions = mmin(num_solutions, N);
    
    double* V = new double[NumSolutions * N];
    double* E = new double[NumSolutions];

    Eigensolver solver;
    solver.SolveLargeSymmetric(M, E, V, N, NumSolutions);

    for(unsigned int i = 0; i < NumSolutions; i++)
    {
        pLevel level(new Level(E[i], (V + N * i), configs, N));
        (*levels)[LevelID(sym, i)] = level;
    }

    delete[] E;
    delete[] V;
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
