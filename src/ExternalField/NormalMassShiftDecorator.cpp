#include "NormalMassShiftDecorator.h"
#include "Include.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"

namespace Ambit
{
NormalMassShiftDecorator::NormalMassShiftDecorator(pHFOperator wrapped_hf, bool only_rel_nms, bool nonrel):
    BaseDecorator(wrapped_hf)
{
    scale = 0.0;
    do_nonrel_nms = !only_rel_nms;
    do_rel_nms = !nonrel;
}

void NormalMassShiftDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    BaseDecorator::GetODEFunction(latticepoint, fg, w);

    if(do_nonrel_nms)
    {
        const double alpha = physicalConstant->GetAlpha();
        const double alphasquared = physicalConstant->GetAlphaSquared();

        w[0] += -2.*scale/alpha * gsl_pow_2(1. + 0.5 * (GetDirectPotential().f[latticepoint] + currentEnergy) * alphasquared) * fg.g[latticepoint];
        w[1] += alpha * 0.5 * scale * gsl_pow_2(alpha * (GetDirectPotential().f[latticepoint] + currentEnergy)) * fg.f[latticepoint];
    }
}

void NormalMassShiftDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    BaseDecorator::GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(do_nonrel_nms)
    {
        const double alpha = physicalConstant->GetAlpha();
        const double alphasquared = physicalConstant->GetAlphaSquared();

        w_g[0] += -2.*scale/alpha * gsl_pow_2(1. + 0.5 * (GetDirectPotential().f[latticepoint] + currentEnergy) * alphasquared);
        w_f[1] += alpha * 0.5 * scale * gsl_pow_2(alpha * (GetDirectPotential().f[latticepoint] + currentEnergy));
    }
}

void NormalMassShiftDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    BaseDecorator::GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(do_nonrel_nms)
    {
        const double alpha = physicalConstant->GetAlpha();
        const double alphasquared = physicalConstant->GetAlphaSquared();

        jacobian[0][1] += -2.*scale/alpha * gsl_pow_2(1. + 0.5 * (GetDirectPotential().f[latticepoint] + currentEnergy) * alphasquared);
        jacobian[1][0] += alpha * 0.5 * scale * gsl_pow_2(alpha * (GetDirectPotential().f[latticepoint] + currentEnergy));
    }
}

SpinorFunction NormalMassShiftDecorator::ApplyTo(const SpinorFunction& a) const
{
    // This will add the relativistic part (via extra exchange) as required
    SpinorFunction ta = BaseDecorator::ApplyTo(a);

    // ApplyTo can be done using the standard formula, since here there is no problem with convergence
    if(do_nonrel_nms)
    {
        const double* R = lattice->R();
        const double* R2 = lattice->Rpower(2);

        std::vector<double> second_derivative_f(a.size());
        std::vector<double> second_derivative_g(a.size());
        std::vector<double> third_derivative_f(a.size());
        std::vector<double> third_derivative_g(a.size());

        differentiator->GetDerivative(a.dfdr, second_derivative_f);
        differentiator->GetDerivative(a.dgdr, second_derivative_g);
        differentiator->GetSecondDerivative(a.dfdr, third_derivative_f);
        differentiator->GetSecondDerivative(a.dgdr, third_derivative_g);

        double coeff_f = a.Kappa() * (a.Kappa() + 1);
        double coeff_g = a.Kappa() * (a.Kappa() - 1);

        for(unsigned int i = 0; i < a.size(); i++)
        {
            ta.f[i] -= 0.5 * scale * (second_derivative_f[i] - coeff_f * a.f[i]/R2[i]);
            ta.dfdr[i] -= 0.5 * scale * (third_derivative_f[i] - coeff_f * (a.dfdr[i] - 2.*a.f[i]/R[i]) /R2[i]);
            ta.g[i] -= 0.5 * scale * (second_derivative_g[i] - coeff_g * a.g[i]/R2[i]);
            ta.dgdr[i] -= 0.5 * scale * (third_derivative_g[i] - coeff_g * (a.dgdr[i] - 2.*a.g[i]/R[i]) /R2[i]);
        }
    }

    return ta;
}

SpinorFunction NormalMassShiftDecorator::CalculateExtraExchange(const SpinorFunction& s) const
{
    SpinorFunction exchange(s.Kappa(), s.size());

    if(do_rel_nms)
    {
        const double* R = lattice->R();
        const double alpha = physicalConstant->GetAlpha();

        double ZalphaOnTwoM = 0.5 * Z * alpha;

        for(int i = 0; i < s.size(); i++)
        {
            exchange.f[i] -= ZalphaOnTwoM/R[i] * (2. * s.dgdr[i] - (s.Kappa()+1) * s.g[i]/R[i]);
            exchange.g[i] += ZalphaOnTwoM/R[i] * (2. * s.dfdr[i] + (s.Kappa()-1) * s.f[i]/R[i]);
        }

        differentiator->GetDerivative(exchange.f, exchange.dfdr);
        differentiator->GetDerivative(exchange.g, exchange.dgdr);
    }

    return exchange;
}

NormalMassShiftOperator::NormalMassShiftOperator(pHFOperator hf, bool only_rel_nms, bool nonrel):
    SpinorOperator(0, hf->GetIntegrator()), hf(hf), do_rel_nms(only_rel_nms), do_nonrel_nms(nonrel)
{
    if(only_rel_nms)
        do_nonrel_nms = false;
    else if(nonrel)
        do_rel_nms = false;
    else
        do_nonrel_nms = do_rel_nms = true;
}

SpinorFunction NormalMassShiftOperator::ReducedApplyTo(const SpinorFunction& a, int kappa_b) const
{
    if(a.Kappa() != kappa_b)
        return SpinorFunction(kappa_b);

    SpinorFunction ret(a.Kappa(), a.size());

    std::vector<double> second_derivative_f(a.size());
    std::vector<double> second_derivative_g(a.size());

    const double* R = hf->GetLattice()->R();
    const double* R2 = hf->GetLattice()->Rpower(2);
    const double alpha = hf->GetPhysicalConstant()->GetAlpha();

    hf->GetDifferentiator()->GetDerivative(a.dfdr, second_derivative_f);
    hf->GetDifferentiator()->GetDerivative(a.dgdr, second_derivative_g);

    if(do_nonrel_nms)
    {
        std::vector<double> third_derivative_f(a.size());
        std::vector<double> third_derivative_g(a.size());

        hf->GetDifferentiator()->GetSecondDerivative(a.dfdr, third_derivative_f);
        hf->GetDifferentiator()->GetSecondDerivative(a.dgdr, third_derivative_g);

        double coeff_f = a.Kappa() * (a.Kappa() + 1);
        double coeff_g = a.Kappa() * (a.Kappa() - 1);

        for(unsigned int i = 0; i < a.size(); i++)
        {
            ret.f[i] -= 0.5 * (second_derivative_f[i] - coeff_f * a.f[i]/R2[i]);
            ret.dfdr[i] -= 0.5 * (third_derivative_f[i] - coeff_f * (a.dfdr[i] - 2.*a.f[i]/R[i]) /R2[i]);
            ret.g[i] -= 0.5 * (second_derivative_g[i] - coeff_g * a.g[i]/R2[i]);
            ret.dgdr[i] -= 0.5 * (third_derivative_g[i] - coeff_g * (a.dgdr[i] - 2.*a.g[i]/R[i]) /R2[i]);
        }
    }

    if(do_rel_nms)
    {
        double ZalphaOnTwoM = 0.5 * hf->GetZ() * alpha;

        for(unsigned int i = 0; i < a.size(); i++)
        {
            ret.f[i] += ZalphaOnTwoM/R[i] * (2. * a.dgdr[i] - (a.Kappa()+1) * a.g[i]/R[i]);
            ret.dfdr[i] += ZalphaOnTwoM/R[i]
                    * (2. * second_derivative_g[i] - (a.Kappa()+3) * a.dgdr[i]/R[i] + (2.*a.Kappa()+2) * a.g[i]/R2[i]);

            ret.g[i] -= ZalphaOnTwoM/R[i] * (2. * a.dfdr[i] + (a.Kappa() - 1) * a.f[i]/R[i]);
            ret.dgdr[i] -= ZalphaOnTwoM/R[i]
                    * (2. * second_derivative_f[i] + (a.Kappa()-3) * a.dfdr[i]/R[i] - (2.*a.Kappa()-2) * a.f[i]/R2[i]);
        }
    }

    // Convert to reduced matrix element
    return ret * std::sqrt(a.TwoJ() + 1);
}

NormalMassShiftCalculator::NormalMassShiftCalculator(MultirunOptions& user_input, Atom& atom):
    TransitionCalculator(user_input, atom.GetBasis(), atom.GetLevels())
{
    auto hf = atom.GetHFOperator();
    op = std::make_shared<NormalMassShiftOperator>(hf, user_input.search("--only-relativistic-nms"), user_input.search("--nonrelativistic-mass-shift"));

    if(user_input.search("--rpa"))
    {   op = MakeRPA(std::static_pointer_cast<NormalMassShiftOperator>(op), hf, atom.GetHartreeY());
        double scale = user_input("Scale", 1.);
        std::static_pointer_cast<RPAOperator>(op)->SetScale(scale);
    }
}

void NormalMassShiftCalculator::PrintHeader() const
{
    if(user_input.search("--reduced-elements"))
        *outstream << "NMS reduced matrix elements (a.u.): " << std::endl;
    else
        *outstream << "NMS matrix elements (stretched states) in a.u.: " << std::endl;
}

void NormalMassShiftCalculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    MathConstant* math = MathConstant::Instance();
    double value = matrix_element;
    if(user_input.search("--reduced-elements"))
    {   int twoj1 = left.first->GetTwoJ();
        int twoj2 = right.first->GetTwoJ();
        value = value/math->Electron3j(twoj2, twoj1, op->GetK(), twoj2, -twoj1);
    }

    *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(8) << value << std::endl;
}

}
