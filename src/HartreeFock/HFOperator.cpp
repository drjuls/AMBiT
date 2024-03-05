#include "HFOperator.h"
#include "Include.h"
#include "Universal/MathConstant.h"

namespace Ambit
{
HFOperator::HFOperator(double Z, pCoreConst hf_core, pPhysicalConstant physical_constant, pIntegrator integration_strategy, pCoulombOperator coulomb) :
    HFOperatorBase(Z, hf_core, physical_constant, integration_strategy), coulombSolver(coulomb), currentExchangePotential(-1)
{}

void HFOperator::SetCore(pCoreConst hf_core)
{
    core = hf_core;

    double N = double(core->NumElectrons());
    charge = Z - N;
    directPotential.resize(lattice->size());
    const double* R = lattice->R();

    // Get electron density function
    RadialFunction density;
    density.Clear();

    for(auto cs = core->begin(); cs != core->end(); cs++)
    {
        const Orbital& s = *cs->second;
        double number_electrons = core->GetOccupancy(&s);

        density += s.GetDensity() * number_electrons;
    }
    
    RadialFunction y(density.size());

    if(density.size())
        coulombSolver->GetPotential(density, y, Z - charge);

    unsigned int i = 0;
    while(i < y.size())
    {   directPotential.f[i] = Z/R[i] - y.f[i];
        directPotential.dfdr[i] = -Z/(R[i]*R[i]) - y.dfdr[i];
        i++;
    }
    while(i < directPotential.size())
    {   directPotential.f[i] = charge/R[i];
        directPotential.dfdr[i] = -charge/(R[i]*R[i]);
        i++;
    }
}

RadialFunction HFOperator::GetDirectPotential() const
{
    return directPotential;
}

void HFOperator::Alert()
{
    unsigned int i = directPotential.size();
    const double* R = lattice->R();

    directPotential.resize(lattice->size());

    while(i < directPotential.size())
    {   directPotential.f[i] = charge/R[i];
        directPotential.dfdr[i] = -charge/(R[i]*R[i]);
        i++;
    }

    if(directPotential.size() < currentExchangePotential.size())
        currentExchangePotential.resize(directPotential.size());
}

void HFOperator::SetODEParameters(int kappa, double energy, const SpinorFunction* exchange)
{
    currentKappa = kappa;
    currentEnergy = energy;

    if(exchange)
        currentExchangePotential = *exchange;
    else
    {   currentExchangePotential.Clear();
        currentExchangePotential.resize(directPotential.size());
    }
}

void HFOperator::SetODEParameters(const Orbital& approximation)
{
    currentExchangePotential = CalculateExchange(approximation);
    currentEnergy = approximation.Energy();
    currentKappa = approximation.Kappa();
}

SpinorFunction HFOperator::GetExchange(pOrbitalConst approximation) const
{
    if(approximation == NULL)
        return currentExchangePotential;

    return CalculateExchange(*approximation);
}

// Get df/dr = w[0] and dg/dr = w[1] given point r, (f, g).
void HFOperator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    const double R = lattice->R(latticepoint);
    const double alpha = physicalConstant->GetAlpha();

    const unsigned int i = latticepoint;

    w[0] = -currentKappa/R * fg.f[i] + (2./alpha + alpha * (currentEnergy + directPotential.f[i])) * fg.g[i];
    w[1] = -alpha * (currentEnergy + directPotential.f[i]) * fg.f[i] + currentKappa/R * fg.g[i];

    if(include_nonlocal && (i < currentExchangePotential.size()))
    {   w[0] = w[0] + alpha * currentExchangePotential.g[i];
        w[1] = w[1] - alpha * currentExchangePotential.f[i];
    }
}

void HFOperator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    const double R = lattice->R(latticepoint);
    const double alpha = physicalConstant->GetAlpha();
    
    const unsigned int i = latticepoint;
    
    w_f[0] = -currentKappa/R;
    w_g[0] = 2./alpha + alpha * (currentEnergy + directPotential.f[i]);
    w_f[1] = -alpha * (currentEnergy + directPotential.f[i]);
    w_g[1] = currentKappa/R;

    if(include_nonlocal && (i < currentExchangePotential.size()))
    {   w_const[0] = alpha * currentExchangePotential.g[i];
        w_const[1] = -alpha * currentExchangePotential.f[i];
    }
    else
    {   w_const[0] = 0.;
        w_const[1] = 0.;
    }
}

// Get Jacobian (dw[i]/df and dw[i]/dg), and dw[i]/dr at a point r, (f, g).
void HFOperator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    const double R = lattice->R(latticepoint);
    const double alpha = physicalConstant->GetAlpha();
    
    const unsigned int i = latticepoint;

    jacobian[0][0] = -currentKappa/R;
    jacobian[0][1] = 2./alpha + alpha * (currentEnergy + directPotential.f[i]);
    jacobian[1][0] = -alpha * (currentEnergy + directPotential.f[i]);
    jacobian[1][1] = currentKappa/R;

    dwdr[0] = currentKappa/(R*R) * fg.f[i] + alpha * directPotential.dfdr[i] * fg.g[i];
    dwdr[1] = -alpha * directPotential.dfdr[i] * fg.f[i] -currentKappa/(R*R) * fg.g[i];

    if(include_nonlocal && (i < currentExchangePotential.size()))
    {   dwdr[0] = dwdr[0] + alpha * currentExchangePotential.dgdr[i];
        dwdr[1] = dwdr[1] -alpha * currentExchangePotential.dfdr[i];
    }
}

void HFOperator::EstimateOrbitalNearOrigin(unsigned int numpoints, SpinorFunction& s) const
{
    const RadialFunction& V = directPotential;

    const int start_point = 0;
    const double alpha = physicalConstant->GetAlpha();

    double correction = 0.;
    if(s.size() >= numpoints)
        correction = s.f[start_point];
    else
        s.resize(numpoints);
    
    unsigned int i;
    for(i=start_point; i<start_point+numpoints; i++)
    {   if(s.Kappa() < 0)
        {   s.f[i] = pow(lattice->R(i), -s.Kappa());
            s.g[i] = alpha * s.f[i] * lattice->R(i) * V.f[i] / (2 * s.Kappa() - 1);
            s.dfdr[i] = - s.Kappa() * s.f[i] / lattice->R(i);
            s.dgdr[i] = ( - s.Kappa() + 1.) * s.g[i] / lattice->R(i);
        }
        else
        {   s.g[i] = alpha * pow(lattice->R(i), s.Kappa());
            s.f[i] = s.g[i] * lattice->R(i) * alpha * V.f[i] / (2 * s.Kappa() + 1);
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

void HFOperator::EstimateOrbitalNearInfinity(unsigned int numpoints, Orbital& s) const
{
    const double alpha = physicalConstant->GetAlpha();

    // Get end point
    unsigned int end_point, original_end_point = 0;
    unsigned int i;
    double correction = 0.;

    if(s.size() < numpoints)
    {   s.resize(directPotential.size());
        end_point = s.size() - 1;
    }
    else
    {   end_point = s.size() - 1;
        original_end_point = end_point;
        correction = s.f[original_end_point];
    }

    i = end_point - numpoints + 1;
    double P;
    double rturn = 0.;
    do
    {   while(i < directPotential.size())
        {   P = -2.*(directPotential.f[i] + s.Energy()) + double(s.Kappa()*(s.Kappa() + 1))/pow(lattice->R(i),2.);
            if(P > 0.)
                break;
            i++;
        }

        if(P < 0.)
        {   // Take action for failing to get to semi-classical regime
            if(rturn == 0.)
            {   // Estimate classical turning point
                rturn = -mmax(GetCharge(), 1.)/s.Energy();
                while(rturn < lattice->MaxRealDistance())
                    rturn *= 1.5;
            }
            else
                rturn *= 1.5;

            lattice->resize_to_r(rturn);
        }
    } while(P < 0.);

    end_point = i + numpoints - 1;
    if(end_point > directPotential.size() - 1)
    {   end_point = directPotential.size() - 1;
    }
    s.resize(end_point+1);

    double S = -9.;
    for(i = end_point; i > end_point - numpoints; i--)
    {
        P = -2.*(directPotential.f[i] + s.Energy()) + s.Kappa()*(s.Kappa() + 1.)/pow(lattice->R(i),2.);
        //assert(P>0);
        P = sqrt(P);
        S = S + 0.5 * P * lattice->dR(i);
        
        s.f[i] = exp(S)/sqrt(P);
        s.g[i] = alpha * s.f[i] * (s.Kappa()/lattice->R(i) - P) * 0.5;
        s.dfdr[i] = (-P * s.f[i]);
        s.dgdr[i] = s.Kappa()/lattice->R(i) * s.g[i] - alpha * (s.Energy() + directPotential.f[i]) * s.f[i];
        
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

SpinorFunction HFOperator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta(a.Kappa(), a.size());

    const double alpha = physicalConstant->GetAlpha();
    const double alphasquared = physicalConstant->GetAlphaSquared();
    const double* R = lattice->R();
    
    double kappa = a.Kappa();
    SpinorFunction exchangePotential(a.Kappa());
    if(include_nonlocal)
        exchangePotential = CalculateExchange(a);
    exchangePotential.resize(ta.size());

    std::vector<double> d2fdr2(a.size());
    std::vector<double> d2gdr2(a.size());
    differentiator->GetDerivative(a.dfdr, d2fdr2);
    differentiator->GetDerivative(a.dgdr, d2gdr2);

    for(unsigned int i = 0; i < ta.size(); i++)
    {
        ta.f[i] = -directPotential.f[i]*a.f[i] - exchangePotential.f[i]
                  + (-a.dgdr[i] + kappa/R[i]*a.g[i])/alpha;
        ta.g[i] = (a.dfdr[i] + kappa/R[i]*a.f[i])/alpha
                  - (2./alphasquared + directPotential.f[i])*a.g[i] - exchangePotential.g[i];

        ta.dfdr[i] = - directPotential.f[i]*a.dfdr[i] - directPotential.dfdr[i]*a.f[i] - exchangePotential.dfdr[i]
                     + (-d2gdr2[i] + kappa/R[i] * (a.dgdr[i] - a.g[i]/R[i]))/alpha;
        ta.dgdr[i] = (d2fdr2[i] + kappa/R[i] * (a.dfdr[i] - a.f[i]/R[i]))/alpha
                     - (2./alphasquared + directPotential.f[i]) * a.dgdr[i] - directPotential.dfdr[i] * a.g[i] - exchangePotential.dgdr[i];
    }

    return ta;
}

SpinorFunction HFOperator::CalculateExchange(const SpinorFunction& s) const
{
    bool NON_REL_SCALING = true;

    SpinorFunction exchange(s.Kappa());
    exchange.resize(s.size());

    // Find out whether s is in the core
    const Orbital* current_in_core = dynamic_cast<const Orbital*>(&s);
    if(current_in_core && core->GetState(OrbitalInfo(current_in_core)) == nullptr)
        current_in_core = nullptr;

    // Sum over all core states
    for(auto cs = core->begin(); cs != core->end(); cs++)
    {
        pOrbitalConst core_orbital = cs->second;
        double other_occupancy = core->GetOccupancy(OrbitalInfo(core_orbital));

        // Get overlap of wavefunctions
        RadialFunction density = s.GetDensity(*core_orbital);

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

            // Integrate density to get (1/r)Y(ab,r)
            RadialFunction potential(mmax(density.size(), core_orbital->size()));
            coulombSolver->GetPotential(k, density, potential);

            exchange += (*core_orbital) * potential * coefficient;
        }
    }

    return exchange;
}
}
