#ifndef HAMILTONIAN_MATRIX_H
#define HAMILTONIAN_MATRIX_H

#include "RelativisticConfiguration.h"
#include "NonRelConfiguration.h"
#include "CIIntegrals.h"
#include "HartreeFock/HFOperator.h"
#include "Level.h"
#include "Basis/ExcitedStates.h"
//#include "MBPT/Sigma3Calculator.h"
#include "Universal/Enums.h"
#include "Universal/Matrix.h"
#include "ManyBodyOperator.h"
#include "OneElectronOperator.h"

typedef OneElectronOperator<pHFOperatorConst> HFElectronOperator;
typedef boost::shared_ptr<HFElectronOperator> pHFElectronOperator;
typedef boost::shared_ptr<const HFElectronOperator> pHFElectronOperatorConst;

typedef ManyBodyOperator<pHFElectronOperator, pTwoElectronCoulombOperator> TwoBodyHamiltonianOperator;
typedef boost::shared_ptr<TwoBodyHamiltonianOperator> pTwoBodyHamiltonianOperator;
typedef boost::shared_ptr<const TwoBodyHamiltonianOperator> pTwoBodyHamiltonianOperatorConst;

class HamiltonianMatrix
{
public:
    HamiltonianMatrix(pHFElectronOperator hf, const CIIntegrals& coulomb_integrals, pRelativisticConfigListConst relconfigs);
    virtual ~HamiltonianMatrix();

    virtual void GenerateMatrix();
    virtual void WriteToFile(const std::string& filename);
    virtual void PollMatrix();

    /** Solve the matrix that has been generated.
        If gFactors are required, set boolean to true.
        min_percentage is the threshold whereby leading configurations are printed
     */
    virtual void SolveMatrix(const Symmetry& sym, unsigned int num_solutions, pLevelMap levels);

    Matrix* GetMatrix()
    {
        return M;
    }

public:
//    /** Include Sigma3 in the leading configurations. */
//    void IncludeSigma3(Sigma3Calculator* sigma3)
//    {   if(sigma3 && (configs->front().ElectronNumber() >= 3))
//            include_sigma3 = true;
//        sigma3calc = sigma3;
//    }
//
//    bool include_sigma3;
//    Sigma3Calculator* sigma3calc;
//    pConfigListConst leading_configs;
//
//    /** Return value of Sigma3 for matrix element (added in GetProjectionH).
//        Checks to see that first and second are leading configurations.
//      */
//    double GetSigma3(const Projection& first, const Projection& second) const;
//
//    /** Get Sigma3(e1 e2 e3 -> e4 e5 e6) including all permutations.
//        Sigma3 does the pair matching (matching e1 to e4, e5, and then e6).
//     */
//    double Sigma3(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
//                  const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const;
//
//    /** This function does the line permutations, putting the pairs on different levels
//        of the three-body interaction.
//     */
//    inline double Sigma3LinePermutations(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
//                  const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const;

protected:
    pRelativisticConfigListConst configs;
    pTwoBodyHamiltonianOperator H_two_body;

    Matrix* M;      // Hamiltonian Matrix
    unsigned int N; // Matrix M dimension = N x N
};

#endif
