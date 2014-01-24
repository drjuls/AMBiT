#include "MassShiftDecorator.h"
#include "Include.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/Interpolator.h"

MassShiftDecorator::MassShiftDecorator(OneBodyOperator* wrapped_OBO, SpinorODE* wrapped_ODE, pOPIntegrator integration_strategy):
    OneBodyOperatorDecorator(wrapped_OBO, integration_strategy), SpinorODEDecorator(wrapped_ODE), extraExchangePotential(-1), lambda(0.0)
{}

/** Set exchange (nonlocal) potential and energy for ODE routines. */
void MassShiftDecorator::SetODEParameters(const SingleParticleWavefunction& approximation)
{
    wrapped->SetODEParameters(approximation);
    extraExchangePotential = CalculateExtraExchange(approximation);
}

/** Get exchange (nonlocal) potential. */
SpinorFunction MassShiftDecorator::GetExchange(pSingleParticleWavefunctionConst approximation) const
{
    SpinorFunction ret = wrapped->GetExchange(approximation);

    if(approximation == NULL)
        ret += extraExchangePotential;
    else
        ret += CalculateExtraExchange(*approximation);

    return ret;
}

void MassShiftDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    wrapped->GetODEFunction(latticepoint, fg, w);

    if(include_nonlocal && latticepoint < extraExchangePotential.Size())
    {   double alpha = PhysicalConstant::Instance()->GetAlpha();
        w[0] += alpha * extraExchangePotential.g[latticepoint];
        w[1] -= alpha * extraExchangePotential.f[latticepoint];
    }
}
void MassShiftDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(include_nonlocal && latticepoint < extraExchangePotential.Size())
    {   double alpha = PhysicalConstant::Instance()->GetAlpha();
        w_const[0] += alpha * extraExchangePotential.g[latticepoint];
        w_const[1] -= alpha * extraExchangePotential.f[latticepoint];
    }
}
void MassShiftDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(include_nonlocal && latticepoint < extraExchangePotential.Size())
    {   double alpha = PhysicalConstant::Instance()->GetAlpha();
        dwdr[0] += alpha * extraExchangePotential.dgdr[latticepoint];
        dwdr[1] -= alpha * extraExchangePotential.dfdr[latticepoint];
    }
}

SpinorFunction MassShiftDecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = component->ApplyTo(a);
    ta -= CalculateExtraExchange(a);
    
    return ta;
}

SpinorFunction MassShiftDecorator::CalculateExtraExchange(const SpinorFunction& s) const
{
    bool NON_REL_SCALING = true;
    
    SpinorFunction exchange(s.Kappa());
    exchange.ReSize(s.Size());
    
    // Find out whether s is in the core
    const Orbital* ds_current = dynamic_cast<const Orbital*>(&s);
    if(core->GetState(OrbitalInfo(ds_current)) == NULL)
        ds_current = NULL;
    
    // Sum over all core states
    ConstStateIterator cs = core->GetConstStateIterator();
    while(!cs.AtEnd())
    {
        pOrbitalConst other = cs.GetState();
        
        // Sum over all k
        for(unsigned int k = abs((int)other->L() - (int)s.L()); k <= (other->L() + s.L()); k+=2)
        {
            double coefficient = MathConstant::Instance()->Electron3j(s.TwoJ(), other->TwoJ(), k);
            coefficient = (2 * abs(other->Kappa())) * coefficient * coefficient;
            
            // Open shells need to be scaled
            if(core->IsOpenShellState(OrbitalInfo(other)) && (other->Occupancy() != double(2 * abs(other->Kappa()))))
            {
                double ex = 1.;
                if(NON_REL_SCALING)
                {   // Average over non-relativistic configurations
                    if(other->Kappa() == -1)
                    {
                        if((ds_current == NULL) || (OrbitalInfo(ds_current) != OrbitalInfo(other)))
                            ex = other->Occupancy()/double(2 * abs(other->Kappa()));
                        else if(k)
                            ex = (other->Occupancy()-1.)/double(2 * abs(other->Kappa()) - 1);
                    }
                    else
                    {
                        int other_kappa = - other->Kappa() - 1;
                        pOrbitalConst ds = core->GetState(OrbitalInfo(other->GetPQN(), other_kappa));
                        
                        if((!ds_current && s.L() != other->L())
                           || (ds_current && (OrbitalInfo(ds_current) != OrbitalInfo(other)) && (OrbitalInfo(ds_current) != OrbitalInfo(ds))))
                            ex = (other->Occupancy() + ds->Occupancy())/double(2 * (abs(other->Kappa()) + abs(ds->Kappa())));
                        else if(k)
                            ex = (other->Occupancy() + ds->Occupancy() - 1.)/double(2 * (abs(other->Kappa()) + abs(ds->Kappa())) - 1);
                    }
                }
                else
                {   // Average over relativistic configurations
                    if((ds_current == NULL) || (OrbitalInfo(ds_current) != OrbitalInfo(other)))
                        ex = other->Occupancy()/double(2 * (abs(other->Kappa())));
                    else if(k)
                        ex = (other->Occupancy() - 1.)/double(2 * (abs(other->Kappa())) - 1);
                }
                
                coefficient = coefficient * ex;
            }

            if(lambda && (k == 1))
            {
                RadialFunction P;
                double sms = CalculateSMS(s, *other, &P);

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
