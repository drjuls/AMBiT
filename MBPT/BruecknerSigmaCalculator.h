#ifndef BRUECKNER_SIGMA_CALCULATOR_H
#define BRUECKNER_SIGMA_CALCULATOR_H

#include "MBPTCalculator.h"
#include "SigmaPotential.h"
#include "HartreeFock/HartreeY.h"

class BruecknerSigmaCalculator : public MBPTCalculator
{
public:
    BruecknerSigmaCalculator(pOrbitalManagerConst orbitals, pSpinorOperatorConst one_body, pHartreeY two_body);
    virtual ~BruecknerSigmaCalculator() {}

    virtual unsigned int GetStorageSize() override { return 0; }
    virtual void UpdateIntegrals() override { return; }

    /** Create a second-order one-electron MBPT (sigma1) operator. */
    void GetSecondOrderSigma(int kappa, SigmaPotential& sigma);

protected:
    /** Calculate diagrams of second order (shown here 1 through 4).
     ->>------>------>>--  ->>------>------<---  ->>------<------>>--  ->>------<------>---
     a   |  beta |   b     a   |  beta |   n     a   |   m   |   b     a   |   m   | alpha
     |       |             |       |             |       |             |       |
     --<------>------<---  --<------>------>>--  --<------>------<---  -->------<------>>--
     n     alpha     n     n     alpha     b     n     alpha     n    alpha    n       b

     PRE: si.kappa == sf.kappa
     */
    void CalculateCorrelation1and3(int kappa, SigmaPotential& sigma);
    void CalculateCorrelation2(int kappa, SigmaPotential& sigma);
    void CalculateCorrelation4(int kappa, SigmaPotential& sigma);

protected:
    pSpinorOperatorConst hf;
    pHartreeY hartreeY;

    pOrbitalMapConst core;
    pOrbitalMapConst excited;
};

typedef std::shared_ptr<BruecknerSigmaCalculator> pBruecknerSigmaCalculator;

#endif
