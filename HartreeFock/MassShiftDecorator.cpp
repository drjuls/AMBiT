#include "MassShiftDecorator.h"
#include "Include.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/Interpolator.h"

MassShiftDecorator::MassShiftDecorator(pHFOperator wrapped_hf, pOPIntegrator integration_strategy):
    HFOperatorDecorator(wrapped_hf), lambda(0.0)
{   if(integration_strategy != NULL)
        integrator = integration_strategy;
}

/** Set exchange (nonlocal) potential and energy for ODE routines. */
void MassShiftDecorator::SetODEParameters(const SingleParticleWavefunction& approximation)
{
    HFOperatorDecorator::SetODEParameters(approximation);
    currentExchangePotential = CalculateExtraExchange(approximation);
}

/** Get exchange (nonlocal) potential. */
SpinorFunction MassShiftDecorator::GetExchange(pSingleParticleWavefunctionConst approximation) const
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

    if(include_nonlocal && latticepoint < currentExchangePotential.Size())
    {   double alpha = PhysicalConstant::Instance()->GetAlpha();
        w[0] += alpha * currentExchangePotential.g[latticepoint];
        w[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}
void MassShiftDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(include_nonlocal && latticepoint < currentExchangePotential.Size())
    {   double alpha = PhysicalConstant::Instance()->GetAlpha();
        w_const[0] += alpha * currentExchangePotential.g[latticepoint];
        w_const[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}
void MassShiftDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(include_nonlocal && latticepoint < currentExchangePotential.Size())
    {   double alpha = PhysicalConstant::Instance()->GetAlpha();
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
    exchange.ReSize(s.Size());

    // Find out whether s is in the core
    const Orbital* current_in_core = dynamic_cast<const Orbital*>(&s);
    if(core->GetState(OrbitalInfo(current_in_core)) == NULL)
        current_in_core = NULL;
    
    // Sum over all core states
    ConstStateIterator cs = core->GetConstStateIterator();
    while(!cs.AtEnd())
    {
        pOrbitalConst core_orbital = cs.GetState();
        double other_occupancy = core->GetOccupancy(OrbitalInfo(core_orbital));

        // Sum over all k
        for(unsigned int k = abs((int)core_orbital->L() - (int)s.L()); k <= (core_orbital->L() + s.L()); k+=2)
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
                    {   OrbitalInfo pair_info(core_orbital->GetPQN(), - core_orbital->Kappa() - 1);
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

            if(lambda && (k == 1))
            {
                RadialFunction P;
                double sms = CalculateSMS(s, *core_orbital, &P);

                for(unsigned int i=0; i < mmin(exchange.Size(), P.Size()); i++)
                {
                    exchange.f[i] += coefficient * lambda * sms *  P.f[i];
                    exchange.dfdr[i] += coefficient * lambda * sms *  P.dfdr[i];
                }
            }

        }
        cs.Next();
    }
    
    return exchange;
}

double MassShiftDecorator::CalculateSMS(const SpinorFunction& s1, const SpinorFunction& s2, RadialFunction* p) const
{
    if(p)
    {   p->Clear();
        p->ReSize(s2.Size());
    }

    double coeff_fb = 0.;

    // Check angular momenta delta functions
    if(s1.L() == s2.L()+1)
        coeff_fb = -double(s1.L());
    else if(s1.L() == s2.L()-1)
        coeff_fb = double(s2.L());
    else
        return 0.;

    // If p is NULL, then make a RadialFunction p_fb, otherwise just use p itself
    RadialFunction* p_fb;
    if(p)
        p_fb = p;
    else
        p_fb = new RadialFunction(s2.Size());

    double total = 0.0;
    const double* R = lattice->R();
    const double* dR = lattice->dR();

    unsigned int i = 0;
    while(i < mmin(s1.Size(), s2.Size()))
    {   p_fb->f[i] = s2.dfdr[i] + coeff_fb/R[i] * s2.f[i];
        total += s1.f[i] * p_fb->f[i] * dR[i];
        i++;
    }

    if(p == NULL)
        delete p_fb;
    else
    {   // Finish making p
        while(i < s2.Size())
        {   p_fb->f[i] = s2.dfdr[i] + coeff_fb/R[i] * s2.f[i];
            i++;
        }

        // Get derivative of p->f
        Interpolator I(lattice);
        I.GetDerivative(s2.dfdr, p->dfdr, 6);  // First term
        // Second term
        for(i = 0; i < p->Size(); i++)
            p->dfdr[i] += coeff_fb/R[i] * (s2.dfdr[i] - s2.f[i]/R[i]);
    }

    return total;
}
