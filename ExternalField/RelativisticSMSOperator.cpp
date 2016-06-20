#include "RelativisticSMSOperator.h"

/** < b | t | a > for an operator t. */
double RelativisticSMSOperator::GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse) const
{
    if(!integrator)
    {   *errstream << "RelativisticSMSOperator::GetMatrixElement(): no integrator found." << std::endl;
        exit(1);
    }

    SpinorFunction ta = ApplyTo(a, b.Kappa(), reverse);
    return integrator->GetInnerProduct(ta, b);
}

SpinorFunction RelativisticSMSOperator::ApplyOperator(const SpinorFunction& a, int kappa_b) const
{
    if(!integrator)
    {   *errstream << "RelativisticSMSOperator::ApplyOperator(): no integrator found." << std::endl;
        exit(1);
    }

    // Check angular momentum selection rules
    if((a.Kappa() != -kappa_b) && (abs(a.Kappa() - kappa_b) != 1))
        return SpinorFunction(kappa_b);

    // Operator p
    SpinorFunction p(kappa_b, a.size());

    auto eta = [](int kappa_i, int kappa_j)
    {
        int l_i = (kappa_i > 0? kappa_i: -kappa_i-1);
        int l_j = (kappa_j > 0? kappa_j: -kappa_j-1);

        switch(l_i - l_j)
        {
            case 1: return -l_i;
            case -1: return l_j;
            default: return 0;
        }
    };

    double eta_fa = eta(kappa_b, a.Kappa());
    double eta_ga = eta(-kappa_b, -a.Kappa());

    pLattice lattice = integrator->GetLattice();
    const double* R = lattice->R();

    unsigned int i = 0;
    for(i = 0; i < a.size(); i++)
    {
        p.f[i] = a.dfdr[i] + eta_fa/R[i] * a.f[i];
        p.g[i] = a.dgdr[i] + eta_ga/R[i] * a.g[i];
    }

    // Get derivative: First term
    interp.GetDerivative(a.dfdr, p.dfdr);
    interp.GetDerivative(a.dgdr, p.dgdr);
    // Second term
    for(i = 0; i < p.size(); i++)
    {   p.dfdr[i] += eta_fa/R[i] * (a.dfdr[i] - a.f[i]/R[i]);
        p.dgdr[i] += eta_ga/R[i] * (a.dgdr[i] - a.g[i]/R[i]);
    }

    if(include_rel)
    {
        // Relativistic part
        auto zeta = [](int kappa_i, int kappa_j)
        {
            if(kappa_i == -kappa_j)
                return 2 * kappa_j + 1;
            else if(kappa_i == kappa_j - 1)
                return 2;
            else
                return 0;
        };

        double zeta_fa = zeta(kappa_b, a.Kappa()) + 1;
        double zeta_ga = zeta(-kappa_b, -a.Kappa()) + 1;

        SpinorFunction r(kappa_b, a.size());

        for(unsigned int i = 0; i < a.size(); i++)
        {
            r.f[i] = - zeta_fa/R[i] * a.g[i];
            r.dfdr[i] = - zeta_fa/R[i] * (a.dgdr[i] - a.g[i]/R[i]);
            r.g[i] = zeta_ga/R[i] * a.f[i];
            r.dgdr[i] = zeta_ga/R[i] * (a.dfdr[i] - a.f[i]/R[i]);
        }

        p += r * (0.5 * Zalpha);
    }

    return p;
}
