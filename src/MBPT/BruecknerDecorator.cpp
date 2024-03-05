#include "BruecknerDecorator.h"

namespace Ambit
{
BruecknerDecorator::BruecknerDecorator(pHFOperator wrapped_hf, pIntegrator integration_strategy):
    HFOperatorDecorator(wrapped_hf, integration_strategy)
{
    SetMatrixParameters();
    IncludeLower();
}

void BruecknerDecorator::SetMatrixParameters(int stride, double start, double end)
{
    matrix_stride = stride;
    matrix_start = lattice->real_to_lattice(start);
    matrix_end = lattice->real_to_lattice(end);
}

void BruecknerDecorator::CalculateSigma(int kappa, pBruecknerSigmaCalculator brueckner_calculator)
{
    if(sigmas.find(kappa) == sigmas.end())
    {
        pSigmaPotential sigma(new SigmaPotential(lattice, matrix_end, matrix_start, matrix_stride));
        sigma->IncludeLower(use_fg, use_gg);
        brueckner_calculator->GetSecondOrderSigma(kappa, *sigma);

        sigmas[kappa] = sigma;
    }
}

void BruecknerDecorator::CalculateSigma(int kappa, pOrbitalManagerConst orbitals, pHartreeY hartreeY, const std::string& fermi_orbitals, pSpinorOperatorConst bare_hf)
{
    if(sigmas.find(kappa) == sigmas.end())
    {
        pBruecknerSigmaCalculator calculator;
        if(bare_hf)
            calculator.reset(new BruecknerSigmaCalculator(orbitals, bare_hf, hartreeY, fermi_orbitals));
        else
            calculator.reset(new BruecknerSigmaCalculator(orbitals, wrapped, hartreeY, fermi_orbitals));

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

SpinorFunction BruecknerDecorator::CalculateExtraExchange(const SpinorFunction& s) const
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

        differentiator->GetDerivative(ret.f, ret.dfdr);
        differentiator->GetDerivative(ret.g, ret.dgdr);
    }

    return ret;
}
}
