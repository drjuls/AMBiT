#include "MassShiftDecorator.h"
#include "NonRelativisticSMSOperator.h"
#include "RelativisticSMSOperator.h"
#include "Include.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"

MassShiftDecorator::MassShiftDecorator(pHFOperator wrapped_hf, bool include_sms, bool include_nms, bool nonrel, bool only_rel_nms, bool include_lower_sms):
    BaseDecorator(wrapped_hf), lambda(0.0), sms_operator(nullptr)
{
    do_nonrel_nms = include_nms && !only_rel_nms;
    do_rel_nms = include_nms && !nonrel;

    if(include_sms)
    {
        pHartreeY zero = std::make_shared<HartreeYBase>();
        if(nonrel && !include_lower_sms)
        {
            sms_operator = std::make_shared<NonRelativisticSMSOperator>(zero, integrator);
        }
        else
        {
            double Zalpha = Z * physicalConstant->GetAlpha();
            sms_operator = std::make_shared<RelativisticSMSOperator>(zero, Zalpha, !nonrel, integrator);
        }

        sms_operator->SetInverseMass(lambda);
    }
}

void MassShiftDecorator::Alert()
{
    if(currentExchangePotential.size() > lattice->size())
        currentExchangePotential.resize(lattice->size());
}

/** Set exchange (nonlocal) potential and energy for ODE routines. */
void MassShiftDecorator::SetODEParameters(const Orbital& approximation)
{
    HFOperatorDecorator::SetODEParameters(approximation);
    currentExchangePotential = CalculateExtraExchange(approximation);
}

/** Get exchange (nonlocal) potential. */
SpinorFunction MassShiftDecorator::GetExchange(pOrbitalConst approximation) const
{
    SpinorFunction ret = wrapped->GetExchange(approximation);

    if(approximation == NULL)
        ret += currentExchangePotential;
    else
        ret += CalculateExtraExchange(*approximation);

    return ret;
}

void MassShiftDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    wrapped->GetODEFunction(latticepoint, fg, w);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        w[0] += alpha * currentExchangePotential.g[latticepoint];
        w[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}
void MassShiftDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        w_const[0] += alpha * currentExchangePotential.g[latticepoint];
        w_const[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}
void MassShiftDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        dwdr[0] += alpha * currentExchangePotential.dgdr[latticepoint];
        dwdr[1] -= alpha * currentExchangePotential.dfdr[latticepoint];
    }
}

SpinorFunction MassShiftDecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = wrapped->ApplyTo(a);
    ta -= CalculateExtraExchange(a);

    return ta;
}

SpinorFunction MassShiftDecorator::CalculateExtraExchange(const SpinorFunction& s) const
{
    bool NON_REL_SCALING = true;

    SpinorFunction exchange(s.Kappa());

    if(lambda == 0)
        return exchange;

    if(sms_operator)
    {
        // Find out whether s is in the core
        const Orbital* current_in_core = dynamic_cast<const Orbital*>(&s);
        if(core->GetState(OrbitalInfo(current_in_core)) == nullptr)
            current_in_core = NULL;

        pSpinorFunctionConst p_s(s.shared_from_this());

        // Sum over all core states
        auto cs = core->begin();
        while(cs != core->end())
        {
            pOrbitalConst core_orbital = cs->second;
            double other_occupancy = core->GetOccupancy(OrbitalInfo(core_orbital));

            int k = sms_operator->SetOrbitals(p_s, core_orbital);

            // Sum over all k
            while(k != -1)
            {
                double coefficient = MathConstant::Instance()->Electron3j(s.TwoJ(), core_orbital->TwoJ(), k);
                coefficient = (2 * abs(core_orbital->Kappa())) * coefficient * coefficient;

                // Open shells need to be scaled
                if(other_occupancy != double(2 * abs(core_orbital->Kappa())))
                {
                    double ex = 1.;
                    if(NON_REL_SCALING)
                    {   // Average over non-relativistic configurations
                        if(core_orbital->Kappa() == -1)
                        {
                            if(!current_in_core || (OrbitalInfo(current_in_core) != OrbitalInfo(core_orbital)))
                                ex = other_occupancy/double(2 * abs(core_orbital->Kappa()));
                            else if(k)
                                ex = (other_occupancy - 1.)/double(2 * abs(core_orbital->Kappa()) - 1);
                        }
                        else
                        {   OrbitalInfo pair_info(core_orbital->PQN(), - core_orbital->Kappa() - 1);
                            pOrbitalConst pair_orbital = core->GetState(pair_info);
                            double pair_occupancy = core->GetOccupancy(pair_info);

                            if((!current_in_core && s.L() != core_orbital->L())
                               || (current_in_core && (OrbitalInfo(current_in_core) != OrbitalInfo(core_orbital)) && (OrbitalInfo(current_in_core) != pair_info)))
                                ex = (other_occupancy + pair_occupancy)/double(2 * (abs(core_orbital->Kappa()) + abs(pair_info.Kappa())));
                            else if(k)
                                ex = (other_occupancy + pair_occupancy - 1.)/double(2 * (abs(core_orbital->Kappa()) + abs(pair_info.Kappa())) - 1);
                        }
                    }
                    else
                    {   // Average over relativistic configurations
                        if(!current_in_core || (OrbitalInfo(current_in_core) != OrbitalInfo(core_orbital)))
                            ex = other_occupancy/double(2 * (abs(core_orbital->Kappa())));
                        else if(k)
                            ex = (other_occupancy - 1.)/double(2 * (abs(core_orbital->Kappa())) - 1);
                    }

                    coefficient = coefficient * ex;
                }

                exchange -= sms_operator->ApplyTo(*core_orbital, s.Kappa()) * coefficient;

                k = sms_operator->NextK();
            }
            cs++;
        }
    }

    if(do_nonrel_nms || do_rel_nms)
    {
        if(exchange.size() < s.size())
            exchange.resize(s.size());

        std::vector<double> second_derivative_f(s.size());
        std::vector<double> second_derivative_g(s.size());

        const double* R = lattice->R();
        const double* R2 = lattice->Rpower(2);

        // Non-relativistic NMS
        differentiator->GetDerivative(s.dfdr, second_derivative_f);
        differentiator->GetDerivative(s.dgdr, second_derivative_g);

        if(do_nonrel_nms)
        {
            std::vector<double> third_derivative_f(s.size());
            std::vector<double> third_derivative_g(s.size());

            differentiator->GetSecondDerivative(s.dfdr, third_derivative_f);
            differentiator->GetSecondDerivative(s.dgdr, third_derivative_g);

            double coeff_f = s.Kappa() * (s.Kappa() + 1);
            double coeff_g = s.Kappa() * (s.Kappa() - 1);

            for(int i = 0; i < s.size(); i++)
            {
                exchange.f[i] += 0.5 * lambda * (second_derivative_f[i] - coeff_f * s.f[i]/R2[i]);
                exchange.dfdr[i] += 0.5 * lambda * (third_derivative_f[i] - coeff_f * (s.dfdr[i] - 2.*s.f[i]/R[i]) /R2[i]);
                exchange.g[i] += 0.5 * lambda * (second_derivative_g[i] - coeff_g * s.g[i]/R2[i]);
                exchange.dgdr[i] += 0.5 * lambda * (third_derivative_g[i] - coeff_g * (s.dgdr[i] - 2.*s.g[i]/R[i]) /R2[i]);
            }
        }

        if(do_rel_nms)
        {
            double ZalphaOnTwoM = 0.5 * lambda * Z * physicalConstant->GetAlpha();

            for(int i = 0; i < s.size(); i++)
            {
                exchange.f[i] -= ZalphaOnTwoM/R[i] * (2. * s.dgdr[i] - (s.Kappa()+1) * s.g[i]/R[i]);
                exchange.dfdr[i] -= ZalphaOnTwoM/R[i]
                    * (2. * second_derivative_g[i] - (s.Kappa()+3) * s.dgdr[i]/R[i] + (2.*s.Kappa()+2) * s.g[i]/R2[i]);

                exchange.g[i] += ZalphaOnTwoM * (2. * s.dfdr[i] + (s.Kappa() - 1) * s.f[i]/R[i])/R[i];
                exchange.dgdr[i] += ZalphaOnTwoM/R[i]
                    * (2. * second_derivative_f[i] + (s.Kappa()-3) * s.dfdr[i]/R[i] - (2.*s.Kappa()-2) * s.f[i]/R2[i]);
            }
        }
    }

    return exchange;
}
