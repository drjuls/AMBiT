#ifndef BRUECKNER_DECORATOR_H
#define BRUECKNER_DECORATOR_H

#include "HartreeFock/HFOperator.h"
#include "SigmaPotential.h"
#include "CoreMBPTCalculator.h"

/** BruecknerDecorator adds one-body "sigma" operators to a HF potential.
    It stores a mapping from kappa to SigmaPotential, and can read and write sigmas to file.
    Creation is deferred to a CoreMBPTCalculator that must be supplied as necessary.

    The sigma potential is non-local, so it is added with the exchange part.
 */
class BruecknerDecorator : public HFOperatorDecorator<BruecknerDecorator>
{
public:
    BruecknerDecorator(pHFOperator wrapped_hf, pOPIntegrator integration_strategy = nullptr);
    BruecknerDecorator(const BruecknerDecorator& other):
        HFOperatorDecorator(other), sigmas(other.sigmas), lambda(other.lambda),
        use_fg(other.use_fg), use_gg(other.use_gg), matrix_stride(other.matrix_stride),
        matrix_start(other.matrix_start), matrix_end(other.matrix_end)
    {}

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
    void CalculateSigma(int kappa, pOrbitalManagerConst orbitals, pHartreeY hartreeY, pSpinorOperatorConst bare_hf = nullptr);

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

public:
    virtual void Alert() override;

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(const Orbital& approximation) override;

    /** Get exchange (nonlocal) potential. */
    virtual SpinorFunction GetExchange(pOrbitalConst approximation) const override;

    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

protected:
    virtual SpinorFunction CalculateExtraNonlocal(const SpinorFunction& s, bool include_derivative) const;

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

#endif
