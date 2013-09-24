#include "HFOperator.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"
#include "Interpolator.h"

HFOperator::HFOperator(double Z, const Core* hf_core, OPIntegrator* integration_strategy) :
    OneBodyOperator(integration_strategy), SpinorODE(hf_core->GetLattice()), currentExchangePotential(-1)
{
    this->Z = Z;
    SetCore(hf_core);
}

HFOperator::~HFOperator()
{}

/** Set/reset the Hartree-Fock core, from which the potential is derived. */
void HFOperator::SetCore(const Core* hf_core)
{
    core = hf_core;
    directPotential = core->GetHFPotential();

    Interpolator interp(core->GetLattice());
    dVdR.resize(directPotential.size());
    interp.GetDerivative(directPotential, dVdR, 6);
}

/** Set exchange (nonlocal) potential and energy for ODE routines. */
void HFOperator::SetODEParameters(int kappa, double energy, SpinorFunction* exchange)
{
    currentKappa = kappa;
    currentEnergy = energy;

    if(exchange)
        currentExchangePotential = *exchange;
    else
    {   currentExchangePotential.Clear();
        currentExchangePotential.ReSize(directPotential.size());
    }
}

void HFOperator::SetODEParameters(const SingleParticleWavefunction* approximation)
{
    core->CalculateExchange(*approximation, currentExchangePotential);
    currentEnergy = approximation->GetEnergy();
    currentKappa = approximation->Kappa();
}

/** Get exchange (nonlocal) potential. */
SpinorFunction HFOperator::GetExchange(const SingleParticleWavefunction* approximation)
{
    if(approximation == NULL)
        return currentExchangePotential;

    SpinorFunction exchange(approximation->Kappa());
    core->CalculateExchange(*approximation, exchange);
    return exchange;
}

/** Get df/dr = w[0] and dg/dr = w[1] given point r, (f, g). */
void HFOperator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    const double R = core->GetLattice()->R(latticepoint);
    const double alpha = PhysicalConstant::Instance()->GetAlpha();

    const unsigned int i = latticepoint;

    w[0] = -currentKappa/R * fg.f[i] + (2./alpha + alpha * (currentEnergy + directPotential[i])) * fg.g[i];
    w[1] = -alpha * (currentEnergy + directPotential[i]) * fg.f[i] + currentKappa/R * fg.g[i];

    if(include_nonlocal)
    {   w[0] = w[0] + alpha * currentExchangePotential.g[i];
        w[1] = w[1] - alpha * currentExchangePotential.f[i];
    }
}

void HFOperator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    const double R = core->GetLattice()->R(latticepoint);
    const double alpha = PhysicalConstant::Instance()->GetAlpha();
    
    const unsigned int i = latticepoint;
    
    w_f[0] = -currentKappa/R;
    w_g[0] = 2./alpha + alpha * (currentEnergy + directPotential[i]);
    w_f[1] = -alpha * (currentEnergy + directPotential[i]);
    w_g[1] = currentKappa/R;

    if(include_nonlocal)
    {   w_const[0] = alpha * currentExchangePotential.g[i];
        w_const[1] = -alpha * currentExchangePotential.f[i];
    }
    else
    {   w_const[0] = 0.;
        w_const[1] = 0.;
    }
}

/** Get Jacobian (dw[i]/df and dw[i]/dg), and dw[i]/dr at a point r, (f, g).
 PRE: jacobian should be an allocated 2x2 matrix,
 dwdr should be an allocated 2 dimensional array. */
void HFOperator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    const double R = core->GetLattice()->R(latticepoint);
    const double alpha = PhysicalConstant::Instance()->GetAlpha();
    
    const unsigned int i = latticepoint;

    jacobian[0][0] = -currentKappa/R;
    jacobian[0][1] = 2./alpha + alpha * (currentEnergy + directPotential[i]);
    jacobian[1][0] = -alpha * (currentEnergy + directPotential[i]);
    jacobian[1][1] = currentKappa/R;

    dwdr[0] = currentKappa/(R*R) * fg.f[i] + alpha * dVdR[i] * fg.g[i];
    dwdr[1] = -alpha * dVdR[i] * fg.f[i] -currentKappa/(R*R) * fg.g[i];

    if(include_nonlocal)
    {   dwdr[0] = dwdr[0] + alpha * currentExchangePotential.dgdr[i];
        dwdr[1] = dwdr[1] -alpha * currentExchangePotential.dfdr[i];
    }
}

/** Get approximation to eigenfunction at a point assumed to be near the origin. */
void HFOperator::EstimateOrbitalNearOrigin(unsigned int numpoints, SpinorFunction& s) const
{
    const int start_point = 0;
    const double alpha = PhysicalConstant::Instance()->GetAlpha();

    double correction = 0.;
    if(s.Size() >= numpoints)
        correction = s.f[start_point];
    else
        s.ReSize(numpoints);
    
    unsigned int i;
    for(i=start_point; i<start_point+numpoints; i++)
    {   if(s.Kappa() < 0)
        {   s.f[i] = pow(lattice->R(i), -s.Kappa());
            s.g[i] = alpha * s.f[i] * lattice->R(i) * directPotential[i] / (2 * s.Kappa() - 1);
            s.dfdr[i] = - s.Kappa() * s.f[i] / lattice->R(i);
            s.dgdr[i] = ( - s.Kappa() + 1.) * s.g[i] / lattice->R(i);
        }
        else
        {   s.g[i] = alpha * pow(lattice->R(i), s.Kappa());
            s.f[i] = s.g[i] * lattice->R(i) * alpha * directPotential[i] / (2 * s.Kappa() + 1);
            s.dgdr[i] = s.Kappa() * s.g[i] / lattice->R(i);
            s.dfdr[i] = (s.Kappa() + 1.) * s.f[i] / lattice->R(i);
        }
    }
    
    // Determine an appropriate scaling to make the norm close to unit.
    if(correction)
        correction = correction/s.f[start_point];
    else
        correction = Z * Z;
    
    for(i=start_point; i<start_point+numpoints; i++)
    {   s.f[i] = s.f[i] * correction;
        s.g[i] = s.g[i] * correction;
        s.dfdr[i] = s.dfdr[i] * correction;
        s.dgdr[i] = s.dgdr[i] * correction;
    }
}

/** Get approximation to eigenfunction at a point assumed to be far from the origin. */
void HFOperator::EstimateOrbitalNearInfinity(unsigned int numpoints, Orbital& s) const
{
    const double alpha = PhysicalConstant::Instance()->GetAlpha();

    // Get end point
    unsigned int end_point, original_end_point = 0;
    unsigned int i;
    double correction = 0.;

    if(s.Size() < numpoints)
    {   s.ReSize(directPotential.size());
        end_point = s.Size() - 1;
    }
    else
    {   end_point = s.Size() - 1;
        original_end_point = end_point;
        correction = s.f[original_end_point];
    }

    i = end_point - numpoints + 1;
    double P;
    while(i < directPotential.size())
    {   P = -2.*(directPotential[i] + s.GetEnergy()) + double(s.Kappa()*(s.Kappa() + 1))/pow(lattice->R(i),2.);
        if(P > 0.)
            break;
        i++;
    }
    
    end_point = i + numpoints - 1;
    if(end_point > directPotential.size() - 1)
    {   end_point = directPotential.size() - 1;
    }
    s.ReSize(end_point+1);

    double S = -9.;
    for(i = end_point; i > end_point - numpoints; i--)
    {
        P = -2.*(directPotential[i] + s.GetEnergy()) + s.Kappa()*(s.Kappa() + 1.)/pow(lattice->R(i),2.);
        //assert(P>0);
        P = sqrt(P);
        S = S + 0.5 * P * lattice->dR(i);
        
        s.f[i] = exp(S)/sqrt(P);
        s.g[i] = alpha * s.f[i] * (s.Kappa()/lattice->R(i) - P) * 0.5;
        s.dfdr[i] = (-P * s.f[i]);
        s.dgdr[i] = s.Kappa()/lattice->R(i) * s.g[i] - alpha * (s.GetEnergy() + directPotential[i]) * s.f[i];
        
        S = S + 0.5 * P * lattice->dR(i);
    }

    if(correction)
    {   correction = correction/s.f[end_point];
        for(i = end_point; i > end_point - numpoints; i--)
        {   s.f[i] *= correction;
            s.g[i] *= correction;
            s.dfdr[i] *= correction;
            s.dgdr[i] *= correction;
        }
    }
}

/** SpinorFunction = t | a > for an operator t. */
SpinorFunction HFOperator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta(a.Size());
    
    const double alpha = PhysicalConstant::Instance()->GetAlpha();
    const double alphasquared = PhysicalConstant::Instance()->GetAlphaSquared();
    const double* R = core->GetLattice()->R();
    
    double kappa = (double) a.Kappa();
    SpinorFunction exchangePotential(a.Kappa());
    core->CalculateExchange(a, exchangePotential);
    
    for(unsigned int i = 0; i < a.Size(); i++)
    {
        ta.f[i] = -directPotential[i]*a.f[i] - exchangePotential.f[i]
        + (-a.dgdr[i] + kappa/R[i]*a.g[i])/alpha;
        ta.g[i] = (a.dfdr[i] + kappa/R[i]*a.f[i])/alpha
        - (2./alphasquared + directPotential[i])*a.g[i] - exchangePotential.g[i];
    }
    
    return ta;
}
