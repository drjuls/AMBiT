#include "BreitZero.h"
#include "Include.h"
#include "Universal/MathConstant.h"

bool BreitZero::SetParameters(int new_K, const SpinorFunction& c, const SpinorFunction& d)
{
    bool component_is_zero = component->SetParameters(new_K, c, d);
    MathConstant& math = *MathConstant::Instance();

    K = new_K;
    pot_P_Kplus.Clear();
    pot_P_Kminus.Clear();
    pot_Q_Kplus.Clear();
    pot_Q_Kminus.Clear();
    pot_V_K.Clear();

    if(K == 0 || !math.triangular_condition(c.TwoJ(), d.TwoJ(), 2*new_K))
    {
        return component_is_zero;
    }

    // M_K(ijkl) and O_K(ijkl)
    double coeff_plus = math.SphericalTensorReducedMatrixElement(c.Kappa(), d.Kappa(), K);
    if(coeff_plus)
    {
        RadialFunction Qcd = Q(c, d);
        coulomb->GetPotential(K+1, Qcd, pot_Q_Kplus);
        pot_Q_Kplus *= coeff_plus;

        RadialFunction Pcd = P(c, d);
        coulomb->GetPotential(K-1, Pcd, pot_P_Kminus);
        pot_P_Kminus *= coeff_plus;

        // Additional for O_K(ijkl)
        coulomb->GetPotential(K-1, Qcd, pot_Q_Kminus);
        coulomb->GetPotential(K+1, Pcd, pot_P_Kplus);
        pot_Q_Kminus *= coeff_plus;
        pot_P_Kplus *= coeff_plus;
    }

    // N_K(ijkl)
    double coeff_minus = math.SphericalTensorReducedMatrixElement(-c.Kappa(), d.Kappa(), K);
    if(coeff_minus)
    {
        RadialFunction Vcd(mmin(c.size(), d.size()));
        for(unsigned int i = 0; i < Vcd.size(); i++)
        {
            Vcd.f[i] = c.f[i] * d.g[i] + c.g[i] * d.f[i];
            Vcd.dfdr[i] = c.f[i] * d.dgdr[i] + c.dfdr[i] * d.g[i] + c.g[i] * d.dfdr[i] + c.dgdr[i] * d.f[i];
        }

        coulomb->GetPotential(K, Vcd, pot_V_K);
        pot_V_K *= coeff_minus * (c.Kappa() + d.Kappa());
    }

    jc_plus_jd_even = (((c.TwoJ() + d.TwoJ())/2)%2 == 0);
    return (pot_Q_Kplus.size() || pot_V_K.size() || !component_is_zero);
}

bool BreitZero::isZero() const
{   return (!pot_Q_Kplus.size() && !pot_V_K.size() && component->isZero());
}

double BreitZero::GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse) const
{
    double value = component->GetMatrixElement(b, a, reverse);
    if(K == 0 || isZero())
        return value;

    MathConstant& math = *MathConstant::Instance();
    double M = 0, N = 0, O = 0;

    double coeff_plus = math.SphericalTensorReducedMatrixElement(b.Kappa(), a.Kappa(), K);
    if(coeff_plus)
    {
        RadialFunction Qba, Pba;

        if(pot_Q_Kplus.size() || pot_P_Kminus.size() || pot_P_Kplus.size())
        {
            Qba = Q(b, a);
        }
        if(pot_P_Kminus.size() || pot_Q_Kplus.size() || pot_Q_Kminus.size())
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

            M *= math.minus_one_to_the_power(K) * coeff_plus;

            // O_K(bcad)
            O = double((K+1)*(K+1))/double((2 * K + 1) * (2 * K + 3)) * QQ_plus;

            double QP_minus = integrator->GetInnerProduct(Qba, pot_P_Kminus);
            double QP_plus  = integrator->GetInnerProduct(Qba, pot_P_Kplus);

            double PQ_minus = integrator->GetInnerProduct(Pba, pot_Q_Kminus);
            double PQ_plus  = integrator->GetInnerProduct(Pba, pot_Q_Kplus);

            O += double(K*K)/double((2 * K + 1) * (2 * K - 1)) * PP_minus;
            O += double(K * (K+1))/double(4 * K + 2) * (QP_minus - QP_plus + PQ_minus - PQ_plus);

            O *= math.minus_one_to_the_power(K+1) * coeff_plus;

            if(reverse && !jc_plus_jd_even)
            {   M = -M;
                O = -O;
            }
        }
    }

    double coeff_minus = math.SphericalTensorReducedMatrixElement(-b.Kappa(), a.Kappa(), K);
    if(coeff_minus && pot_V_K.size())
    {
        RadialFunction Vba(mmin(b.size(), a.size()));
        for(unsigned int i = 0; i < Vba.size(); i++)
        {
            Vba.f[i] = b.f[i] * a.g[i] + b.g[i] * a.f[i];
            Vba.dfdr[i] = b.f[i] * a.dgdr[i] + b.dfdr[i] * a.g[i] + b.g[i] * a.dfdr[i] + b.dgdr[i] * a.f[i];
        }

        double VV = integrator->GetInnerProduct(Vba, pot_V_K);
        N = math.minus_one_to_the_power(K+1) * double(b.Kappa() + a.Kappa())/double(K * (K+1)) * VV * coeff_minus;

        if(reverse && jc_plus_jd_even)
            N = -N;
    }

    return (value + M + N + O);
}

SpinorFunction BreitZero::ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse) const
{
    SpinorFunction ret = component->ApplyTo(a, kappa_b, reverse);
    if(K == 0 || isZero())
        return ret;

    MathConstant& math = *MathConstant::Instance();

    double coeff_plus = math.SphericalTensorReducedMatrixElement(kappa_b, a.Kappa(), K);
    if(coeff_plus)
    {
        SpinorFunction Qa(kappa_b), Pa(kappa_b);

        if(pot_Q_Kplus.size() || pot_P_Kminus.size() || pot_P_Kplus.size())
        {
            Qa = Q(a, kappa_b);
        }
        if(pot_P_Kminus.size() || pot_Q_Kplus.size() || pot_Q_Kminus.size())
        {
            Pa = P(a, kappa_b);
        }

        if(Qa.size() || Pa.size())
        {
            SpinorFunction QQ_plus = Qa * pot_Q_Kplus;
            SpinorFunction PP_minus = Pa * pot_P_Kminus;

            double M_coeff = math.minus_one_to_the_power(K) * coeff_plus;
            double O_coeff = math.minus_one_to_the_power(K+1) * coeff_plus;

            if(reverse && !jc_plus_jd_even)
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
            ret += Qa * pot_P_Kminus * QP_coeff;
            ret += Qa * pot_P_Kplus * -QP_coeff;
            ret += Pa * pot_Q_Kminus * QP_coeff;
            ret += Pa * pot_Q_Kplus * -QP_coeff;
        }
    }

    double coeff_minus = math.SphericalTensorReducedMatrixElement(-kappa_b, a.Kappa(), K);
    if(coeff_minus && pot_V_K.size())
    {
        double N_coeff = math.minus_one_to_the_power(K+1) * double(kappa_b + a.Kappa())/double(K * (K+1)) * coeff_minus;
        if(reverse && jc_plus_jd_even)
            N_coeff = -N_coeff;

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
