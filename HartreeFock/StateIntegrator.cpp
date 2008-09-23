#include "Include.h"
#include "StateIntegrator.h"
#include "Universal/Constant.h"

void StateIntegrator::IntegrateForwards(State& s, const std::vector<double>& HFPotential, const CoupledFunction* exchange, int end_point, double nuclear_charge)
{
    const int start_point = 0;

    // Set initial conditions
    StateFunction A(lattice);
    A.SetHFPotential(HFPotential);
    A.SetState(&s);
    A.SetExchange(exchange);

    SetUpForwardsIntegral(s, HFPotential, nuclear_charge);

    Integrate2(A, s, start_point+(adams_N-1), end_point);
}

void StateIntegrator::IntegrateBackwards(State& s, const std::vector<double>& HFPotential, const CoupledFunction* exchange, int end_point)
{
    // Set initial conditions
    StateFunction A(lattice);
    A.SetHFPotential(HFPotential);
    A.SetState(&s);
    A.SetExchange(exchange);

    SetUpBackwardsIntegral(s, HFPotential);

    Integrate2(A, s, s.Size()-adams_N, end_point);
}

unsigned int StateIntegrator::IntegrateBackwardsUntilPeak(State& s, const std::vector<double>& HFPotential, int end_point)
{
    // Set initial conditions
    StateFunction A(lattice);
    A.SetHFPotential(HFPotential);
    A.SetState(&s);

    SetUpBackwardsIntegral(s, HFPotential);

    int i = s.Size()-(adams_N-1);    // Starting point

    // Do one step at a time until peak is reached
    while((s.df[i]/s.df[i+1] > 0.) && (i>end_point))
    {   i--;
        Integrate2(A, s, i, i-1);
    }

    return i;
}

void StateIntegrator::SetUpForwardsIntegral(State& s, const std::vector<double>& HFPotential, double nuclear_charge)
{
    const int start_point = 0;

    double correction = s.f[start_point];

    unsigned int i;
    for(i=start_point; i<start_point+(adams_N-1); i++)
    {   if(s.Kappa() < 0)
        {   s.f[i] = pow(lattice->R(i), -s.Kappa());
            s.g[i] = s.f[i] * lattice->R(i) * HFPotential[i] / (2 * s.Kappa() - 1);
            s.df[i] = - s.Kappa() * s.f[i] / lattice->R(i) * lattice->dR(i);
            s.dg[i] = ( - s.Kappa() + 1.) * s.g[i] / lattice->R(i) * lattice->dR(i);
        }
        else
        {   s.g[i] = pow(lattice->R(i), s.Kappa());
            s.f[i] = s.g[i] * lattice->R(i) * Constant::AlphaSquared * HFPotential[i] / (2 * s.Kappa() + 1);
            s.dg[i] = s.Kappa() * s.g[i] / lattice->R(i) * lattice->dR(i);
            s.df[i] = (s.Kappa() + 1.) * s.f[i] / lattice->R(i) * lattice->dR(i);
        }
    }

    // Determine an appropriate scaling to make the norm close to unit.
    if(s.f[start_point] && correction)
        correction = correction/s.f[start_point];
    else
        correction = nuclear_charge * nuclear_charge;

    for(i=start_point; i<start_point+(adams_N-1); i++)
    {   s.f[i] = s.f[i] * correction;
        s.g[i] = s.g[i] * correction;
        s.df[i] = s.df[i] * correction;
        s.dg[i] = s.dg[i] * correction;
    }
}

void StateIntegrator::SetUpBackwardsIntegral(State& s, const std::vector<double>& HFPotential)
{
    // Get start point
    unsigned int start_point = s.Size() - 1;
    unsigned int i = start_point - (adams_N-2);
    double P;
    while(i < HFPotential.size())
    {   P = -2.*(HFPotential[i] + s.Energy()) + double(s.Kappa()*(s.Kappa() + 1))/pow(lattice->R(i),2.);
        if(P > 0.)
            break;
        i++;
    }

    start_point = i + (adams_N-2);
    if(start_point > HFPotential.size() - 1)
    {   start_point = HFPotential.size() - 1;
    }
    s.ReSize(start_point+1);

    double correction = s.f[start_point];
    double S = -9.;
    for(i=start_point; i>start_point-(adams_N-1); i--)
    {
        P = -2*(HFPotential[i] + s.Energy()) + s.Kappa()*(s.Kappa() + 1)/pow(lattice->R(i),2.);
        //assert(P>0);
        P = sqrt(P);
        S = S + 0.5 * P * lattice->dR(i);

        s.f[i] = exp(S)/sqrt(P);
        s.g[i] = s.f[i] * (s.Kappa()/lattice->R(i) - P) * 0.5;
        s.df[i] = (-P * s.f[i]) * lattice->dR(i);
        s.dg[i] = (s.Kappa()/lattice->R(i) * s.g[i] - (s.Energy() + HFPotential[i]) * s.f[i]) * lattice->dR(i);

        S = S + 0.5 * P * lattice->dR(i);
    }

    if(correction)
    {   correction = correction/s.f[start_point];
        for(unsigned int i=start_point; i>start_point-(adams_N-1); i--)
        {   s.f[i] = s.f[i] * correction;
            s.g[i] = s.g[i] * correction;
            s.df[i] = s.df[i] * correction;
            s.dg[i] = s.dg[i] * correction;
        }
    }
}

unsigned int StateIntegrator::IntegrateContinuum(ContinuumState& s, const std::vector<double>& HFPotential, const CoupledFunction& exchange, double nuclear_charge, double accuracy, double& final_amplitude, double& final_phase)
{
    // Start calculation outside of the nucleus
    unsigned int start_point = lattice->real_to_lattice(1.e-4);

    StateFunction A(lattice);
    A.SetHFPotential(HFPotential);
    A.SetState(&s);
    A.SetExchange(&exchange);

    double correction = s.f[start_point];

    // Set up inital conditions
    SetUpContinuum(s, HFPotential, A, nuclear_charge, start_point);

    if(correction)
    {   correction = correction/s.f[start_point];
        for(unsigned int i = start_point; i < start_point+(adams_N-1); i++)
        {   s.f[i] = s.f[i] * correction;
            s.g[i] = s.g[i] * correction;
            s.df[i] = s.df[i] * correction;
            s.dg[i] = s.dg[i] * correction;
        }
    }

    // Do integration
    std::vector<double> P(HFPotential.size());
    std::vector<double> S(HFPotential.size());
    const double* R = lattice->R();
    const double* dR = lattice->dR();
    unsigned int start_sine = 0;
    final_amplitude = 0.;
    double peak_phase = 0.;

    int i = start_point+(adams_N-1);
    while(i < s.Size())
    {
        double Pot = HFPotential[i] - double(s.Kappa()*(s.Kappa()+1))/(2*R[i]*R[i]);
        P[i] = sqrt(2. * fabs(s.Energy() + Pot));
        S[i] = S[i-1];
        for(unsigned int j=0; j<adams_N; j++)
            S[i] = S[i] + adams_coeff[j] * P[i-j] * dR[i-j];

        if(!start_sine)
        {
            Integrate2(A, s, i, i+1);
            
            if((HFPotential[i] < 2.5 * s.Energy()) && (s.df[i-1]/s.df[i-2] < 0.))
            {
                double x_max;
                double f_max = FindExtremum(s.f, i-2, x_max);
                bool maximum = (s.df[i-2] > 0.);
                double p_max = FindExtremum(P, i-2, x_max, false);

                double old_amplitude = final_amplitude;
                final_amplitude = fabs(f_max * sqrt(p_max));
                if(fabs((final_amplitude - old_amplitude)/final_amplitude) < accuracy)
                {
                    double old_peak_phase = peak_phase;
                    peak_phase = FindExtremum(S, i-2, x_max, false);
                    if((old_peak_phase != 0.) && (fabs((peak_phase - old_peak_phase)/Constant::Pi - 1.) < accuracy))
                    {
                        start_sine = i;
                        if(!maximum)
                            peak_phase = peak_phase - Constant::Pi;
                    }
                }
                else
                    peak_phase = 0.;    // Reset
            }
        }
        else
        {   // quasiclassical region - sine wave approximation
            double amplitude = final_amplitude/sqrt(P[i]);
            double phase = S[i] - peak_phase;
            s.f[i] = amplitude * cos(phase);
            s.g[i] = 0.5 * amplitude * (-P[i] * sin(phase) + s.Kappa()/R[i] * cos(phase));
        }
        i++;
    }

    if(start_sine)
    {
        GetDerivative(s.f, s.df, start_sine, s.Size()-2);
        GetDerivative(s.g, s.dg, start_sine, s.Size()-2);
        GetDerivativeEnd(s.f, s.df, s.Size());
        GetDerivativeEnd(s.g, s.dg, s.Size());

        final_phase = S[s.Size() - 1] - peak_phase;
    }

    return start_sine;
}

double StateIntegrator::HamiltonianMatrixElement(const State& s1, const State& s2, const Core& core)
{
    double E = 0.;
    const double core_pol = core.GetPolarisability();

    if(s1.Kappa() == s2.Kappa())
    {
        CoupledFunction exchange;
        core.CalculateExchange(s2, exchange);

        std::vector<double> Potential(core.GetHFPotential());

        const double* R = lattice->R();
        const double* dR = lattice->dR();

        unsigned int i;
        for(i=0; i<mmin(s1.Size(), s2.Size()); i++)
        {
            if(core_pol)
            {   double r4 = R[i]*R[i] + core.GetClosedShellRadius()*core.GetClosedShellRadius();
                r4 *= r4;

                Potential[i] -= 0.5 * core.GetPolarisability()/r4;
            }

            double EQ1 = s2.df[i]/dR[i] + s2.Kappa()*s2.f[i]/R[i]
                        - (2. + Constant::AlphaSquared*Potential[i])*s2.g[i]
                        - Constant::AlphaSquared * exchange.g[i];
            double EQ2 = s2.dg[i]/dR[i] - s2.Kappa()*s2.g[i]/R[i] + Potential[i]*s2.f[i] + exchange.f[i];

            E = E - (s1.f[i] * EQ2 - s1.g[i] * EQ1)*dR[i];
//            *errstream << i << "\t" << EQ1/s2.g[i]/Constant::AlphaSquared
//                            << "\t" << -EQ2/s2.f[i]
//                            << "\t" << s1.g[i] * EQ1 * dR[i]
//                            << "\t" << -s1.f[i] * EQ2 * dR[i]
//                            << "\t" << -(s1.f[i] * EQ2 - s1.g[i] * EQ1)*dR[i] << std::endl;
        }
    }
    return E;
}

void StateIntegrator::SetUpContinuum(ContinuumState& s, const std::vector<double>& HFPotential, const StateFunction& state_function, double nuclear_charge, unsigned int start_point)
{
    double& Z = nuclear_charge;
    double energy = s.Energy() - Z/lattice->R(start_point) + HFPotential[start_point];
    double AM = 1./sqrt(2. * fabs(energy));
    double GAM = sqrt(s.Kappa()*s.Kappa() - Constant::AlphaSquared*Z*Z);
    double GAM1 = 2.*GAM + 1.;

    double E, ALAMBD;
    if(energy < 0.)
    {   ALAMBD = sqrt(1. - 0.25*Constant::AlphaSquared/(AM*AM))/AM;
        E = 1. - 0.5*Constant::AlphaSquared/(AM*AM);
    }
    else
    {   ALAMBD = sqrt(1. + 0.25*Constant::AlphaSquared/(AM*AM))/AM;
        E = 1. + 0.5*Constant::AlphaSquared/(AM*AM);
    }

    double AQ = pow(2.*Z, GAM)/sqrt(Z * pow(s.Nu(), 3.));
    AQ = -AQ * abs(s.Kappa())/s.Kappa();

    if(energy < 0.)
    {
        double RN = Z*E/ALAMBD - GAM;
        double AQA = (GAM - s.Kappa())/(AM*ALAMBD*(GAM - s.Kappa() + Z*(1. - E)/ALAMBD));
        double AQQ = AQ * AQA;
        for(unsigned int i = start_point; i<start_point+(adams_N-1); i++)
        {
            double Y = 2.*ALAMBD * lattice->R(i);
            double Fre1 = GipReal(-RN, GAM1, Y);
            double Fre2 = GipReal(1.-RN, GAM1, Y);

            double AQI = AQQ * exp(-0.5 * Y) * pow(lattice->R(i), GAM);
            s.f[i] = AQI * AM * ALAMBD * ((Z/ALAMBD - s.Kappa()) * Fre1 - RN * Fre2);
            s.g[i] = -AQI * 0.5/AM * ((Z/ALAMBD - s.Kappa()) * Fre1 + RN * Fre2);

            s.df[i] = (state_function.Coeff1(i) * s.f[i] + state_function.Coeff2(i) * s.g[i] + state_function.Coeff3(i)) * lattice->dR(i);
            s.dg[i] = (state_function.Coeff4(i) * s.f[i] + state_function.Coeff5(i) * s.g[i] + state_function.Coeff6(i)) * lattice->dR(i);
        }
    }
    else
    {
        double CSI = 0.5 * atan2(Z/ALAMBD * (E*s.Kappa() - GAM), (Z/ALAMBD)*(Z/ALAMBD)*E + GAM*s.Kappa());
        double ANA = Z*E/ALAMBD;
        for(unsigned int i = start_point; i<start_point+(adams_N-1); i++)
        {
            double Y = 2.*ALAMBD * lattice->R(i);
            std::complex<double> A(GAM, -ANA), B(GAM1, 0.), C(0., -Y);
            std::complex<double> F = Gip(A, B, C);

            std::complex<double> CS = std::polar(1., 0.5 * Y + CSI);
            std::complex<double> FG = CS * F;

            double AQQ = pow(lattice->R(i), GAM) * 2. * Z * AM * AQ;
            s.f[i] = -AM * ALAMBD * AQQ * FG.imag();
            s.g[i] = -0.5 / AM * AQQ * FG.real();

            s.df[i] = (state_function.Coeff1(i) * s.f[i] + state_function.Coeff2(i) * s.g[i] + state_function.Coeff3(i)) * lattice->dR(i);
            s.dg[i] = (state_function.Coeff4(i) * s.f[i] + state_function.Coeff5(i) * s.g[i] + state_function.Coeff6(i)) * lattice->dR(i);
        }
    }
}

double StateIntegrator::FindExtremum(const std::vector<double>& function, unsigned int zero_point, double& x, bool calculate_x)
{
    double f1 = function[zero_point-1];
    double f2 = function[zero_point];
    double f3 = function[zero_point+1];
    double f4 = function[zero_point+2];

    double A = (f4 - 3*f3 + 3*f2 - f1)/(6 * pow(lattice->H(), 3.));
    double B = (f3 - 2*f2 + f1)/(2 * lattice->H() * lattice->H());
    double C = (-f4 + 6*f3 - 3*f2 - 2*f1)/(6 * lattice->H());
    double D = f2;

    if(calculate_x)
    {   x = 0.;
        if(fabs(B/D) > 1.e-10)
        {   if(fabs(A/D) > 1.e-10)
                x = -B/(3.*A)*(1.-sqrt(1. - 3.*A*C/(B*B)));
            else
                x = -C/(2*B);
        }
    }
    return A*x*x*x + B*x*x + C*x + D;
}

std::complex<double> StateIntegrator::Gip(std::complex<double> A, std::complex<double> B, std::complex<double> Z)
{
    std::complex<double> F(1., 0.);
    std::complex<double> R(1., 0.);

    B = std::conj(B);
    for(unsigned int i=0; i < 10000; i++)
    {
        std::complex<double> X = R * A * B * Z;
        R = X /(std::abs(B)*std::abs(B) * (i+1));
        F = F + R;

        if(std::abs(R)/std::abs(F) < 0.0001)
            break;

        A = A + 1.;
        B = B + 1.;
    }

    return F;
}

double StateIntegrator::GipReal(double A, double B, double Z)
{
    double F = 1.;
    double R = 1.;

    for(unsigned int i=0; i< 1000; i++)
    {
        double X = R * A * B * Z;
        R = X/(B*B * (i+1));
        F = F + R;

        if(fabs(R/F) < 0.001)
            break;

        A++;
        B++;
    }

    return F;
}

double StateIntegrator::StateFunction::Coeff1(int point) const
{
    return  -kappa/lattice->R(point);
}

double StateIntegrator::StateFunction::Coeff2(int point) const
{
    return (2. + Constant::AlphaSquared * (energy + (*HFPotential)[point]));
}

double StateIntegrator::StateFunction::Coeff3(int point) const
{
    if(exchange && (unsigned int)point < exchange->Size())
        return Constant::AlphaSquared * exchange->g[point];
    else return 0.;
}

double StateIntegrator::StateFunction::Coeff4(int point) const
{
    return -(energy + (*HFPotential)[point]);
}

double StateIntegrator::StateFunction::Coeff5(int point) const
{
    return kappa/lattice->R(point);
}

double StateIntegrator::StateFunction::Coeff6(int point) const
{
    if(exchange && (unsigned int)point < exchange->Size())
        return -exchange->f[point];
    else return 0.;
}

double StateIntegrator::IsotopeShiftIntegral(const std::vector<double> f, unsigned int L, const State& s2, std::vector<double>* P)
{
    double coeff_f2;
    if(L == s2.L() + 1)
    {   coeff_f2 = - double(L);
    }
    else
    {   coeff_f2 = double(s2.L());
    }

    const double* R = lattice->R();
    const double* dR = lattice->dR();

    double SMS = 0.;
    if(P == NULL)
        for(unsigned int i=0; i<mmin(f.size(), s2.Size()); i++)
        {
            SMS += f[i] * (s2.df[i] + coeff_f2 * s2.f[i]/R[i] * dR[i]);
        }
    else
    {   unsigned int state_limit = mmin(f.size(), s2.Size());
        unsigned int P_limit = mmin(s2.Size(), P->size());

        unsigned int i;
        for(i=0; i<mmin(P_limit, state_limit); i++)
        {   (*P)[i] = s2.df[i]/dR[i] + coeff_f2 * s2.f[i]/R[i];
            SMS += f[i] * (*P)[i] * dR[i];
        }
        while(i < state_limit)
        {   SMS += f[i] * (s2.df[i] + coeff_f2 * s2.f[i]/R[i] * dR[i]);
            i++;
        }
        while(i < P_limit)
        {   (*P)[i] = s2.df[i]/dR[i] + coeff_f2 * s2.f[i]/R[i];
            i++;
        }
    }

    return SMS;
}

double StateIntegrator::IsotopeShiftIntegral(const State& s1, const State& s2, std::vector<double>* P)
{
    return IsotopeShiftIntegral(s1.f, s1.L(), s2, P);
}

void StateIntegrator::IsotopeShiftIntegral(unsigned int L, const State& s2, std::vector<double>* P)
{
    double coeff_f2;
    if(L == s2.L() + 1)
    {   coeff_f2 = - double(L);
    }
    else
    {   coeff_f2 = double(s2.L());
    }

    const double* R = lattice->R();
    const double* dR = lattice->dR();

    unsigned int P_limit = mmin(s2.Size(), P->size());
    for(unsigned int i = 0; i < P_limit; i++)
        (*P)[i] = s2.df[i]/dR[i] + coeff_f2 * s2.f[i]/R[i];
}
