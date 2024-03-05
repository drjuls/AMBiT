#ifndef BRUECKNER_DECORATOR_H
#define BRUECKNER_DECORATOR_H

#include "HartreeFock/ExchangeDecorator.h"
#include "SigmaPotential.h"
#include "BruecknerSigmaCalculator.h"

namespace Ambit
{
/** BruecknerDecorator adds one-body "sigma" operators to a HF potential.
    It stores a mapping from kappa to SigmaPotential, and can read and write sigmas to file.
    Creation is deferred to a CoreMBPTCalculator that must be supplied as necessary.

    The sigma potential is non-local, so it is added with the exchange part.
 */
class BruecknerDecorator : public HFOperatorDecorator<ExchangeDecorator, BruecknerDecorator>
{
public:
    BruecknerDecorator(pHFOperator wrapped_hf, pIntegrator integration_strategy = nullptr);

    /** Set parameters for new sigma matrices: matrix size = (start_i - end_i)/stride
        where start = lattice->R[start_i] etc.
        Since many outer products are required to make the matrix, smaller matrices are much more efficient.
     */
    void SetMatrixParameters(int stride = 4, double start = 4.35e-5, double end = 8.0);

    void IncludeLower(bool include_fg = false, bool include_gg = false)
    {   use_fg = include_fg;
        use_gg = include_gg;
    }

    /** Calculate sigma for given kappa if it does not already exist. */
    void CalculateSigma(int kappa, pBruecknerSigmaCalculator brueckner_calculator);

    /** Calculate sigma for given kappa if it does not already exist.
        If bare_hf is not given, then use decorated HFOperator.
     */
    void CalculateSigma(int kappa, pOrbitalManagerConst orbitals, pHartreeY hartreeY, const std::string& fermi_orbitals = "", pSpinorOperatorConst bare_hf = nullptr);
    inline void CalculateSigma(int kappa, pOrbitalManagerConst orbitals, pHartreeY hartreeY, pSpinorOperatorConst bare_hf)
    {   CalculateSigma(kappa, orbitals, hartreeY, "", bare_hf);
    }

    /** Set the scaling parameter for given kappa. */
    void SetSigmaScaling(int kappa, double sigma_amount) { lambda[kappa] = sigma_amount; }

    /** Get the scaling parameter for given kappa. */
    double GetSigmaScaling(int kappa) const;

    /** Attempt to read sigma with given kappa, filename is "identifier.[kappa].sigma". */
    void Read(const std::string& identifier, int kappa);

    /** Write sigmas, filename "identifier.[kappa].sigma". */
    void Write(const std::string& identifier, int kappa) const;

    /** Write all sigmas; filenames "identifier.[kappa].sigma". */
    void Write(const std::string& identifier) const;

protected:
    virtual SpinorFunction CalculateExtraExchange(const SpinorFunction& s) const override;

protected:
    std::map<int, pSigmaPotential> sigmas;  //!< Map kappa to Sigma
    std::map<int, double> lambda;           //!< Map kappa to scalings
    bool use_fg, use_gg;
    int matrix_stride;
    int matrix_start;
    int matrix_end;
};

typedef std::shared_ptr<BruecknerDecorator> pBruecknerDecorator;
typedef std::shared_ptr<const BruecknerDecorator> pBruecknerDecoratorConst;

}
#endif
