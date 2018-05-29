#include "TwoBodySMSOperator.h"
#include "Include.h"

namespace Ambit
{
bool TwoBodySMSOperator::SetLocalParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d)
{
    K = new_K;
    c = new_c;
    d = new_d;

    if(K == 1 && lambda && (c->L() + d->L())%2)
    {
        if(!lightweight_mode)
        {
            SpinorFunction td = ApplyOperator(*new_d, new_c->Kappa());
            p_cd = integrator->GetInnerProduct(td, *new_c) * (-lambda);
        }
        else
            p_cd = 1.0;
    }
    else
    {   p_cd = 0.0;
    }

    return (p_cd != 0.0);
}

bool TwoBodySMSOperator::isZeroLocal() const
{
    return (K != 1 || p_cd == 0);
}

int TwoBodySMSOperator::GetLocalMinK() const
{
    if((c->Kappa() == -d->Kappa()) ||
       (abs(c->Kappa() - d->Kappa()) == 1))
        return 1;
    else
        return -1;
}

int TwoBodySMSOperator::GetLocalMaxK() const
{
    if((c->Kappa() == -d->Kappa()) ||
       (abs(c->Kappa() - d->Kappa()) == 1))
        return 1;
    else
        return -1;
}

double TwoBodySMSOperator::GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse) const
{
    if(!integrator)
    {   *errstream << "TwoBodySMSOperator::GetMatrixElement(): no integrator found." << std::endl;
        exit(1);
    }

    SpinorFunction ta = ApplyTo(a, b.Kappa(), reverse);
    return integrator->GetInnerProduct(ta, b);
}

SpinorFunction TwoBodySMSOperator::ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse) const
{
    if(K == 1 && p_cd && !lightweight_mode)
    {
        SpinorFunction ret = ApplyOperator(a, kappa_b) * p_cd;
        if(reverse)
            ret = ret * (-1.0);

        ret += component->ApplyTo(a, kappa_b, reverse);

        return ret;
    }
    else
        return component->ApplyTo(a, kappa_b, reverse);
}

SpinorFunction TwoBodySMSOperator::ApplyOperator(const SpinorFunction& a, int kappa_b) const
{
    if(!integrator)
    {   *errstream << "TwoBodySMSOperator::ApplyOperator(): no integrator found." << std::endl;
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

    if(include_lower)
    {
        unsigned int i = 0;
        for(i = 0; i < a.size(); i++)
        {   p.f[i] = a.dfdr[i] + eta_fa/R[i] * a.f[i];
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
    }
    else
    {
        unsigned int i = 0;
        for(i = 0; i < a.size(); i++)
        {   p.f[i] = a.dfdr[i] + eta_fa/R[i] * a.f[i];
        }

        // Get derivative: First term
        interp.GetDerivative(a.dfdr, p.dfdr);
        // Second term
        for(i = 0; i < p.size(); i++)
        {   p.dfdr[i] += eta_fa/R[i] * (a.dfdr[i] - a.f[i]/R[i]);
        }
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
}
