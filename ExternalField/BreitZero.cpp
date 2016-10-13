#include "BreitZero.h"
#include "Include.h"
#include "Universal/MathConstant.h"

bool BreitZero::isZeroLocal() const
{
    return (!pot_Q_Kplus.size() && !pot_V_K.size());
}

int BreitZero::GetLocalMinK() const
{
    return mmax(abs(c->TwoJ() - d->TwoJ())/2, 1);
}

int BreitZero::GetLocalMaxK() const
{
    return (c->TwoJ() + d->TwoJ())/2;
}

bool BreitZero::SetLocalParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d)
{
    K = new_K;
    c = new_c;
    d = new_d;

    MathConstant& math = *MathConstant::Instance();

    pot_P_Kminus.Clear();
    pot_Q_Kplus.Clear();
    pot_V_K.Clear();
    pot_P_fwd.Clear();
    pot_Q_back.Clear();

    if(K == 0 || !math.triangular_condition(c->TwoJ(), d->TwoJ(), 2*K))
    {
        return false;
    }

    int pot_size = integrator->GetLattice()->size();

    // M_K(ijkl) and O_K(ijkl)
    if(math.sum_is_even(c->L(), d->L(), K))
    {
        // Q potentials
        RadialFunction Qcd = Q(*c, *d);
        pot_Q_Kplus.resize(pot_size);
        pot_Q_back.resize(pot_size);
        RadialFunction temp(pot_size);
        coulomb->GetBackwardPotential(K+1, Qcd, pot_Q_back);
        coulomb->GetForwardPotential(K+1, Qcd, pot_Q_Kplus);
        pot_Q_Kplus += pot_Q_back;

        coulomb->GetBackwardPotential(K-1, Qcd, temp);
        pot_Q_back = temp - pot_Q_back;

        // P potentials
        RadialFunction Pcd = P(*c, *d);
        pot_P_Kminus.resize(pot_size);
        pot_P_fwd.resize(pot_size);
        coulomb->GetForwardPotential(K-1, Pcd, pot_P_fwd);
        coulomb->GetBackwardPotential(K-1, Pcd, pot_P_Kminus);
        pot_P_Kminus += pot_P_fwd;

        coulomb->GetForwardPotential(K+1, Pcd, temp);
        pot_P_fwd -= temp;
    }
    else // N_K(ijkl)
    {
        RadialFunction Vcd(mmin(c->size(), d->size()));
        for(unsigned int i = 0; i < Vcd.size(); i++)
        {
            Vcd.f[i] = c->f[i] * d->g[i] + c->g[i] * d->f[i];
            Vcd.dfdr[i] = c->f[i] * d->dgdr[i] + c->dfdr[i] * d->g[i] + c->g[i] * d->dfdr[i] + c->dgdr[i] * d->f[i];
        }

        pot_V_K.resize(pot_size);
        coulomb->GetPotential(K, Vcd, pot_V_K);
        pot_V_K *= double(c->Kappa() + d->Kappa());
    }

    return (pot_Q_Kplus.size() || pot_V_K.size());
}

double BreitZero::GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse) const
{
    double value = component->GetMatrixElement(b, a, reverse);

    MathConstant& math = *MathConstant::Instance();
    if(K == 0 || isZero() || !math.triangular_condition(a.TwoJ(), b.TwoJ(), 2*K))
        return value;

    double M = 0., N = 0., O = 0.;

    if(math.sum_is_even(a.L(), b.L(), K))
    {
        RadialFunction Qba, Pba;

        if(pot_Q_Kplus.size() || pot_P_Kminus.size() || pot_P_fwd.size())
        {
            Qba = Q(b, a);
        }
        if(pot_P_Kminus.size() || pot_Q_Kplus.size() || pot_Q_back.size())
        {
            Pba = P(b, a);
        }

        if(Qba.size() || Pba.size())
        {
            double QQ_plus = integrator->GetInnerProduct(Qba, pot_Q_Kplus);
            double PP_minus = integrator->GetInnerProduct(Pba, pot_P_Kminus);

            // M_K(bcad)
            M = double(K+1)/double(2 * K + 3) * QQ_plus;
            M += double(K)/double(2 * K - 1) * PP_minus;

            // O_K(bcad)
            O = - double((K+1)*(K+1))/double((2 * K + 1) * (2 * K + 3)) * QQ_plus;
            O -= double(K*K)/double((2 * K + 1) * (2 * K - 1)) * PP_minus;


            double QP_fwd = integrator->GetInnerProduct(Qba, pot_P_fwd);
            double PQ_back = integrator->GetInnerProduct(Pba, pot_Q_back);

            O -= double(K * (K+1))/double(4 * K + 2) * (QP_fwd + PQ_back);

            if(reverse)
            {   M = -M;
                O = -O;
            }
        }
    }
    else if(pot_V_K.size())
    {
        RadialFunction Vba(mmin(b.size(), a.size()));
        for(unsigned int i = 0; i < Vba.size(); i++)
        {
            Vba.f[i] = b.f[i] * a.g[i] + b.g[i] * a.f[i];
            Vba.dfdr[i] = b.f[i] * a.dgdr[i] + b.dfdr[i] * a.g[i] + b.g[i] * a.dfdr[i] + b.dgdr[i] * a.f[i];
        }

        double VV = integrator->GetInnerProduct(Vba, pot_V_K);
        N = - double(b.Kappa() + a.Kappa())/double(K * (K+1)) * VV;
    }

    return (value + M + N + O);
}

SpinorFunction BreitZero::ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse) const
{
    SpinorFunction ret = component->ApplyTo(a, kappa_b, reverse);

    MathConstant& math = *MathConstant::Instance();
    int b_twoJ = 2 * abs(kappa_b) - 1;
    int b_L = (kappa_b > 0? kappa_b: -kappa_b-1);

    if(K == 0 || isZero() || !math.triangular_condition(a.TwoJ(), b_twoJ, 2*K))
        return ret;

    if(math.sum_is_even(a.L(), b_L, K))
    {
        SpinorFunction Qa(kappa_b), Pa(kappa_b);

        if(pot_Q_Kplus.size() || pot_P_Kminus.size() || pot_P_fwd.size())
        {
            Qa = Q(a, kappa_b);
        }
        if(pot_P_Kminus.size() || pot_Q_Kplus.size() || pot_Q_back.size())
        {
            Pa = P(a, kappa_b);
        }

        if(Qa.size() || Pa.size())
        {
            SpinorFunction QQ_plus = Qa * pot_Q_Kplus;
            SpinorFunction PP_minus = Pa * pot_P_Kminus;

            int M_coeff = 1;
            int O_coeff = -1;

            if(reverse)
            {   M_coeff = -M_coeff;
                O_coeff = -O_coeff;
            }

            // M_K(bcad)
            ret += QQ_plus * (double(K+1)/double(2 * K + 3) * M_coeff);
            ret += PP_minus * (double(K)/double(2 * K - 1) * M_coeff);

            // O_K(bcad)
            ret += QQ_plus * (double((K+1)*(K+1))/double((2 * K + 1) * (2 * K + 3)) * O_coeff);
            ret += PP_minus * (double(K*K)/double((2 * K + 1) * (2 * K - 1)) * O_coeff);

            double QP_coeff = double(K * (K+1))/double(4 * K + 2) * O_coeff;
            ret += Qa * pot_P_fwd * QP_coeff;
            ret += Pa * pot_Q_back * QP_coeff;
        }
    }
    else if(pot_V_K.size())
    {
        double N_coeff = - double(kappa_b + a.Kappa())/double(K * (K+1));

        SpinorFunction Va(kappa_b);
        Va.f = a.g;
        Va.dfdr = a.dgdr;
        Va.g = a.f;
        Va.dgdr = a.dfdr;

        ret += Va * pot_V_K * N_coeff;
    }
    
    return ret;
}

SpinorFunction BreitZero::P(const SpinorFunction& k, int kappa_i) const
{
    SpinorFunction ret(kappa_i);
    ret.resize(k.size());

    double upper_factor = 1. + double(k.Kappa() - kappa_i)/K;
    double lower_factor = -1. + double(k.Kappa() - kappa_i)/K;

    for(unsigned int i = 0; i < k.size(); i++)
    {
        ret.f[i]    = upper_factor * k.g[i];
        ret.dfdr[i] = upper_factor * k.dgdr[i];
        ret.g[i]    = lower_factor * k.f[i];
        ret.dgdr[i] = lower_factor * k.dfdr[i];
    }

    return ret;
}

RadialFunction BreitZero::P(const SpinorFunction& i, const SpinorFunction& k) const
{
    RadialFunction ret = i.GetDensity(P(k, i.Kappa()));
    return ret;
}

SpinorFunction BreitZero::Q(const SpinorFunction& k, int kappa_i) const
{
    SpinorFunction ret(kappa_i);
    ret.resize(k.size());

    double upper_factor = -1. + double(k.Kappa() - kappa_i)/(K+1);
    double lower_factor = 1. + double(k.Kappa() - kappa_i)/(K+1);

    for(unsigned int i = 0; i < k.size(); i++)
    {
        ret.f[i]    = upper_factor * k.g[i];
        ret.dfdr[i] = upper_factor * k.dgdr[i];
        ret.g[i]    = lower_factor * k.f[i];
        ret.dgdr[i] = lower_factor * k.dfdr[i];
    }

    return ret;
}

RadialFunction BreitZero::Q(const SpinorFunction& i, const SpinorFunction& k) const
{
    RadialFunction ret = i.GetDensity(Q(k, i.Kappa()));
    return ret;
}
