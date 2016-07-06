#include "LocalPotentialDecorator.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"
#include "CoulombOperator.h"
#include "Universal/Interpolator.h"

RadialFunction LocalPotentialDecorator::GetDirectPotential() const
{
    RadialFunction ret = wrapped->GetDirectPotential();
    ret += directPotential * scale;

    return ret;
}

void LocalPotentialDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    wrapped->GetODEFunction(latticepoint, fg, w);

    if(latticepoint < directPotential.size())
    {   const double alpha = physicalConstant->GetAlpha() * scale;
        w[0] += alpha * directPotential.f[latticepoint] * fg.g[latticepoint];
        w[1] -= alpha * directPotential.f[latticepoint] * fg.f[latticepoint];
    }
}

void LocalPotentialDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(latticepoint < directPotential.size())
    {   const double alpha = physicalConstant->GetAlpha() * scale;
        w_g[0] += alpha * directPotential.f[latticepoint];
        w_f[1] -= alpha * directPotential.f[latticepoint];
    }
}

void LocalPotentialDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(latticepoint < directPotential.size())
    {   const double alpha = physicalConstant->GetAlpha() * scale;
        jacobian[0][1] += alpha * directPotential.f[latticepoint];
        jacobian[1][0] -= alpha * directPotential.f[latticepoint];

        dwdr[0] += alpha * directPotential.dfdr[latticepoint] * fg.g[latticepoint];
        dwdr[1] -= alpha * directPotential.dfdr[latticepoint] * fg.f[latticepoint];
    }
}

void LocalPotentialDecorator::EstimateOrbitalNearOrigin(unsigned int numpoints, SpinorFunction& s) const
{
    RadialFunction V(wrapped->GetDirectPotential());
    V += directPotential * scale;

    const int start_point = 0;
    const double alpha = physicalConstant->GetAlpha();

    double correction = 0.;
    if(s.size() >= numpoints)
        correction = s.f[start_point];
    else
        s.resize(numpoints);

    unsigned int i;
    for(i=start_point; i<start_point+numpoints; i++)
    {
        if(s.Kappa() < 0)
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

SpinorFunction LocalPotentialDecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = wrapped->ApplyTo(a);
    ta -= a * directPotential * scale;
    
    return ta;
}

ImportedPotentialDecorator::ImportedPotentialDecorator(pHFOperator wrapped_hf, const std::string& filename, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy)
{
    std::ifstream infile(filename.c_str());
    if(!infile.is_open())
    {   *errstream << "ImportedPotentialDecorator: unable to open file " << filename << std::endl;
        return;
    }
    std::vector<double> lat, pot;

    double x, y;
    while(infile.good())
    {
        infile >> x >> y;

        if(!infile.fail())
        {
            lat.push_back(x);
            pot.push_back(y);
        }
    }
    infile.close();

    Interpolator interp(lat);
    const double* R = lattice->R();
    double maximum_distance = mmin(lattice->MaxRealDistance(), x);
    directPotential.Clear();
    directPotential.resize(lattice->real_to_lattice(maximum_distance));

    for(int i = 0; i < directPotential.size(); i++)
        interp.Interpolate(pot, R[i], directPotential.f[i], directPotential.dfdr[i], 6);
}

LocalExchangeApproximation::LocalExchangeApproximation(pHFOperator wrapped_hf, pCoulombOperator coulomb, double x_alpha, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy), coulombSolver(coulomb), Xalpha(x_alpha)
{}

void LocalExchangeApproximation::SetCore(pCoreConst hf_core)
{
    HFOperatorDecorator::SetCore(hf_core);
    double min_charge = !charge;
    const double* R = lattice->R();

    // Get electron density function
    RadialFunction density;

    auto cs = core->begin();
    while(cs != core->end())
    {
        const Orbital& s = *cs->second;
        double number_electrons = core->GetOccupancy(OrbitalInfo(&s));
        
        density += s.GetDensity() * number_electrons;
        
        cs++;
    }

    RadialFunction y(density.size());

    if(density.size())
        coulombSolver->GetPotential(density, y, Z-charge);

    unsigned int i = 0;

    // Get local exchange approximation
    directPotential.resize(lattice->size());
    double C = Xalpha * 0.635348143228;  // (81/32\pi^2)^(1/3)

    for(i = 0; i < density.size(); i++)
    {
        directPotential.f[i] = C * pow((density.f[i]/(R[i]*R[i])), 1./3.);
        directPotential.dfdr[i] = C/3. * pow((density.f[i]/(R[i]*R[i])), -2./3.) * (density.dfdr[i] - 2.*density.f[i]/R[i])/(R[i]*R[i]);

        if(directPotential.f[i] + Z/R[i] - y.f[i] < min_charge/R[i])
            break;
    }

    while(i < density.size())
    {   directPotential.f[i] = (-Z + charge + min_charge)/R[i] + y.f[i];
        directPotential.dfdr[i] = -(-Z + charge + min_charge)/(R[i]*R[i]) + y.dfdr[i];
        i++;
    }

    while(i < lattice->size())
    {   directPotential.f[i] = min_charge/R[i];
        directPotential.dfdr[i] = - min_charge/(R[i]*R[i]);
        i++;
    }
}

void LocalExchangeApproximation::Alert()
{
    if(directPotential.size() > lattice->size())
        directPotential.resize(lattice->size());
}
