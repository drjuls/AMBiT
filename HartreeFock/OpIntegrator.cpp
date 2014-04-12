//
//  OpIntegrator.cpp
//  AMBiT
//
//  Created by Julian Berengut on 16/08/13.
//  Copyright (c) 2013 UNSW. All rights reserved.
//

#include "OpIntegrator.h"

double OPIntegrator::GetInnerProduct(const SpinorFunction& a, const SpinorFunction& b) const
{
    if(a.Kappa() != b.Kappa())
        return 0.;

    RadialFunction integrand = a.GetDensity(b);
    return Integrate(integrand);
}

double SimpsonsIntegrator::Integrate(const RadialFunction& integrand) const
{
    double total = 0.;
    const double* dR = lattice->dR();

    if(integrand.size() == 0)
        return 0.0;

    unsigned int i;
    for(i = 1; i < integrand.size()-1; i+=2)
    {
        total += 4. * integrand.f[i] * dR[i]
                + 2. * integrand.f[i+1] * dR[i+1];
    }
    total = total/3.;
    total += integrand.f[0] * dR[0];

    while(i < integrand.size())
    {   total += integrand.f[i] * dR[i];
        i++;
    }

    return total;
}