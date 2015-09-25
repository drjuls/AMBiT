#include "BruecknerDecorator.h"
#include "Universal/Interpolator.h"

BruecknerDecorator::BruecknerDecorator(pHFOperator wrapped_hf, pOPIntegrator integration_strategy):
    HFOperatorDecorator(wrapped_hf, integration_strategy)
{
    IncludeLower();
}

BruecknerDecorator* BruecknerDecorator::Clone() const
{
    pHFOperator wrapped_clone(wrapped->Clone());
    BruecknerDecorator* ret = new BruecknerDecorator(wrapped_clone, integrator);
    ret->sigmas = sigmas;
    ret->use_fg = use_fg;
    ret->use_gg = use_gg;
    ret->lambda = lambda;
    return ret;
}

void BruecknerDecorator::CalculateSigma(int kappa, pBruecknerSigmaCalculator brueckner_calculator)
{
    if(sigmas.find(kappa) == sigmas.end())
    {
        pSigmaPotential sigma(new SigmaPotential(lattice, lattice->real_to_lattice(8.), 100));
        sigma->IncludeLower(use_fg, use_gg);
        brueckner_calculator->GetSecondOrderSigma(kappa, *sigma);

        sigmas[kappa] = sigma;
    }
}

void BruecknerDecorator::CalculateSigma(int kappa, pOrbitalManagerConst orbitals, pHartreeY hartreeY, pSpinorOperatorConst bare_hf)
{
    if(sigmas.find(kappa) == sigmas.end())
    {
        pBruecknerSigmaCalculator calculator;
        if(bare_hf)
            calculator.reset(new BruecknerSigmaCalculator(orbitals, bare_hf, hartreeY));
        else
            calculator.reset(new BruecknerSigmaCalculator(orbitals, wrapped, hartreeY));

        CalculateSigma(kappa, calculator);
    }
}

double BruecknerDecorator::GetSigmaScaling(int kappa) const
{
    double scaling = 1.;

    auto it = lambda.find(kappa);
    if(it != lambda.end())
        scaling = it->second;

    return scaling;
}

/** Attempt to read sigma with given kappa, filename is "identifier.[kappa].sigma". */
void BruecknerDecorator::Read(const std::string& identifier, int kappa)
{
    std::string filename = identifier + "." + itoa(kappa) + ".sigma";

    pSigmaPotential sigma(new SigmaPotential(lattice));
    if(sigma->Read(filename))
    {
        sigmas[kappa] = sigma;
    }
}

/** Write all sigmas; filenames "identifier.[kappa].sigma". */
void BruecknerDecorator::Write(const std::string& identifier, int kappa) const
{
    auto it = sigmas.find(kappa);
    if(it != sigmas.end())
    {
        std::string filename = identifier + "." + itoa(kappa) + ".sigma";
        it->second->Write(filename);
    }
}

/** Write all sigmas; filenames "identifier.[kappa].sigma". */
void BruecknerDecorator::Write(const std::string& identifier) const
{
    for(auto& pair: sigmas)
    {
        std::string filename = identifier + "." + itoa(pair.first) + ".sigma";
        pair.second->Write(filename);
    }
}

void BruecknerDecorator::Alert()
{
    // Don't even try to change the size of the Sigma matrices.
    if(currentExchangePotential.size() > lattice->size())
        currentExchangePotential.resize(lattice->size());
}

/** Set exchange (nonlocal) potential and energy for ODE routines. */
void BruecknerDecorator::SetODEParameters(const Orbital& approximation)
{
    HFOperatorDecorator::SetODEParameters(approximation);
    currentExchangePotential = CalculateExtraNonlocal(approximation, true);
}

/** Get exchange (nonlocal) potential. */
SpinorFunction BruecknerDecorator::GetExchange(pOrbitalConst approximation) const
{
    SpinorFunction ret = wrapped->GetExchange(approximation);

    if(approximation == NULL)
        ret += currentExchangePotential;
    else
        ret += CalculateExtraNonlocal(*approximation, true);

    return ret;
}

void BruecknerDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    wrapped->GetODEFunction(latticepoint, fg, w);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        w[0] += alpha * currentExchangePotential.g[latticepoint];
        w[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}

void BruecknerDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        w_const[0] += alpha * currentExchangePotential.g[latticepoint];
        w_const[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}

void BruecknerDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        dwdr[0] += alpha * currentExchangePotential.dgdr[latticepoint];
        dwdr[1] -= alpha * currentExchangePotential.dfdr[latticepoint];
    }
}

SpinorFunction BruecknerDecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = wrapped->ApplyTo(a);
    ta -= CalculateExtraNonlocal(a, false);

    return ta;
}

SpinorFunction BruecknerDecorator::CalculateExtraNonlocal(const SpinorFunction& s, bool include_derivative) const
{
    SpinorFunction ret(s.Kappa());

    auto it = sigmas.find(s.Kappa());
    if(it != sigmas.end())
    {
        pSigmaPotential sigma = it->second;
        if(s.size() < sigma->size())
        {   SpinorFunction bigger_s(s);
            bigger_s.resize(sigma->size());
            ret = sigma->ApplyTo(bigger_s);
        }
        else
            ret = sigma->ApplyTo(s);

        double scaling = GetSigmaScaling(s.Kappa());
        ret *= (-scaling);  // Remember we are storing HF potential -V (that is, V > 0)

        if(include_derivative)
        {
            Interpolator interp(lattice);
            interp.GetDerivative(ret.f, ret.dfdr, 6);
            interp.GetDerivative(ret.g, ret.dgdr, 6);
        }
    }

    return ret;
}
