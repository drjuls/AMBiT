#include "Include.h"
#include "Integrator.h"

static const double adams_10[10]
= {2082753./7257600., 9449717./7257600., -11271304./7257600., 16002320./7257600., -17283646./7257600.,
   13510082./7257600., -7394032./7257600., 2687864./7257600., -583435./7257600., 57281./7257600.};

static const double adams_5[5]
= {251./720., 646./720., -264./720., 106./720., -19./720.};


void Integrator::SetAdamsOrder(unsigned int adams_order)
{
    adams_N = adams_order;

    if(adams_N <= 5)
    {   adams_N = 5;
        adams_coeff = adams_5;
    }
    else
    {   adams_N = 10;
        adams_coeff = adams_10;
    }
}

void Integrator::Integrate0(const std::vector<double>& f, std::vector<double>& y,
                       int start_point, int end_point)
{
    const double* dR = lattice->dR();

    if(start_point < end_point)
        for(int i = start_point; i < end_point; i++)
        {
            double y_next = y[i-1]
                            + adams_coeff[0] * f[i] * dR[i];

            for(unsigned int j=1; j<adams_N; j++)
            {
                y_next = y_next + adams_coeff[j] * f[i-j] * dR[i-j];
            }

            y[i] = y_next;
            // dy[i] = f[i] * dR[i];
        }
    else
        for(int i = start_point; i > end_point; i--)
        {
            double y_next = y[i+1]
                            - adams_coeff[0] * f[i] * dR[i];

            for(unsigned int j=1; j<adams_N; j++)
            {
                y_next = y_next - adams_coeff[j] * f[i+j] * dR[i+j];
            }

            y[i] = y_next;
            // dy[i] = f[i] * dR[i];
        }
}

void Integrator::Integrate(const Function2& f, std::vector<double>& y, std::vector<double>& dy,
                       int start_point, int end_point)
{
    std::vector<double> Coeff1 = f.Coeff1();
    std::vector<double> Coeff2 = f.Coeff2();
    const double* dR = lattice->dR();

    if(start_point < end_point)
        for(int i = start_point; i < end_point; i++)
        {
            double y_next = y[i-1]
                            + adams_coeff[0] * Coeff2[i] * dR[i];

            for(unsigned int j=1; j<adams_N; j++)
            {
                y_next = y_next + adams_coeff[j] * dy[i-j];
            }

            y[i] = y_next/(1. - adams_coeff[0] * Coeff1[i] * dR[i]);
            dy[i] = (Coeff1[i] * y[i] + Coeff2[i]) * dR[i];
        }
    else
        for(int i = start_point; i > end_point; i--)
        {
            double y_next = y[i+1]
                            - adams_coeff[0] * Coeff2[i] * dR[i];

            for(unsigned int j=1; j<adams_N; j++)
            {
                y_next = y_next - adams_coeff[j] * dy[i+j];
            }

            y[i] = y_next/(1. + adams_coeff[0] * Coeff1[i] * dR[i]);
            dy[i] = (Coeff1[i] * y[i] + Coeff2[i]) * dR[i];
        }
}

void Integrator::Integrate2(const Function6& A, CoupledFunction& s, int start_point, int end_point)
{
    if(start_point < end_point)
        for(int i = start_point; i < end_point; i++)
        {
            double f_next = s.f[i-1] + adams_coeff[0] * A.Coeff3(i) * lattice->dR(i);
            double g_next = s.g[i-1] + adams_coeff[0] * A.Coeff6(i) * lattice->dR(i);
            
            for(unsigned int j=1; j<adams_N; j++)
            {
                f_next = f_next + adams_coeff[j] * s.df[i-j];
                g_next = g_next + adams_coeff[j] * s.dg[i-j];
            }
            
            double D = (1. - adams_coeff[0] * A.Coeff1(i) * lattice->dR(i))
                         * (1. - adams_coeff[0] * A.Coeff5(i) * lattice->dR(i))
                     - (adams_coeff[0] * A.Coeff2(i) * lattice->dR(i))
                         * (adams_coeff[0] * A.Coeff4(i) * lattice->dR(i));

            s.f[i] = (f_next * (1. - adams_coeff[0] * A.Coeff5(i) * lattice->dR(i))
                         + g_next * (adams_coeff[0] * A.Coeff2(i) * lattice->dR(i))) / D;
            s.g[i] = (g_next * (1. - adams_coeff[0] * A.Coeff1(i) * lattice->dR(i))
                         + f_next * (adams_coeff[0] * A.Coeff4(i) * lattice->dR(i))) / D;
            s.df[i] = (A.Coeff1(i) * s.f[i] + A.Coeff2(i) * s.g[i] + A.Coeff3(i)) * lattice->dR(i);
            s.dg[i] = (A.Coeff4(i) * s.f[i] + A.Coeff5(i) * s.g[i] + A.Coeff6(i)) * lattice->dR(i);
        }
    else
        for(int i = start_point; i > end_point; i--)
        {
            double f_next = s.f[i+1] - adams_coeff[0] * A.Coeff3(i) * lattice->dR(i);
            double g_next = s.g[i+1] - adams_coeff[0] * A.Coeff6(i) * lattice->dR(i);
            
            for(unsigned int j=1; j<adams_N; j++)
            {
                f_next = f_next - adams_coeff[j] * s.df[i+j];
                g_next = g_next - adams_coeff[j] * s.dg[i+j];
            }
            
            double D = (1. + adams_coeff[0] * A.Coeff1(i) * lattice->dR(i))
                         * (1. + adams_coeff[0] * A.Coeff5(i) * lattice->dR(i))
                     - (adams_coeff[0] * A.Coeff2(i) * lattice->dR(i))
                         * (adams_coeff[0] * A.Coeff4(i) * lattice->dR(i));

            s.f[i] = (f_next * (1. + adams_coeff[0] * A.Coeff5(i) * lattice->dR(i))
                         - g_next * (adams_coeff[0] * A.Coeff2(i) * lattice->dR(i))) / D;
            s.g[i] = (g_next * (1. + adams_coeff[0] * A.Coeff1(i) * lattice->dR(i))
                         - f_next * (adams_coeff[0] * A.Coeff4(i) * lattice->dR(i))) / D;
            s.df[i] = (A.Coeff1(i) * s.f[i] + A.Coeff2(i) * s.g[i] + A.Coeff3(i)) * lattice->dR(i);
            s.dg[i] = (A.Coeff4(i) * s.f[i] + A.Coeff5(i) * s.g[i] + A.Coeff6(i)) * lattice->dR(i);
        }
}

void Integrator::GetDerivative(const std::vector<double>& f, std::vector<double>& df, int start_point, int end_point)
{
    const double* R = lattice->R();

    // Calculate derivative using quartic formula
    for(int i = start_point; i<end_point; i++)
    {
        double x = R[i];
        double a = R[i-2];
        double b = R[i-1];
        double c = R[i];
        double d = R[i+1];
        double e = R[i+2];

        double coeff_A = (x-b)*(x-d)*(x-e)/((a-b)*(a-c)*(a-d)*(a-e));
        double coeff_B = (x-a)*(x-d)*(x-e)/((b-a)*(b-c)*(b-d)*(b-e));
        double coeff_C = ((x-a)*(x-b)*(x-d) + (x-a)*(x-b)*(x-e) + (x-a)*(x-d)*(x-e) + (x-b)*(x-d)*(x-e))
                          /((c-a)*(c-b)*(c-d)*(c-e));
        double coeff_D = (x-a)*(x-b)*(x-e)/((d-a)*(d-b)*(d-c)*(d-e));
        double coeff_E = (x-a)*(x-b)*(x-d)/((e-a)*(e-b)*(e-c)*(e-d));

        double A = f[i-2];
        double B = f[i-1];
        double C = f[i];
        double D = f[i+1];
        double E = f[i+2];

        df[i] = coeff_A * A + coeff_B * B + coeff_C * C + coeff_D * D + coeff_E * E;
        df[i] = df[i] * lattice->dR(i);
    }
}

/** Calculate the first two points of the derivative of f:
        df[start_point], df[start_point+1].
    Requires f to be known from (start_point) to (start_point+3)
    */
void Integrator::GetDerivativeStart(const std::vector<double>& f, std::vector<double>& df, int start_point)
{
    const double* R = lattice->R();

    // first two points just use quadratic formula
    for(int i = start_point; i < start_point+2; i++)
    {
        double x = R[i];
        double c = R[i];
        double d = R[i+1];
        double e = R[i+2];

        double coeff_C = ((x-d) + (x-e))/((c-d)*(c-e));
        double coeff_D = (x-e)/((d-c)*(d-e));
        double coeff_E = (x-d)/((e-c)*(e-d));

        double C = f[i];
        double D = f[i+1];
        double E = f[i+2];

        df[i] = coeff_C * C + coeff_D * D + coeff_E * E;
        df[i] = df[i] * lattice->dR(i);
    }
}

/** Calculate the last two points of the derivative of f:
        df[end_point-2], df[end_point-1].
    Requires f to be known from (end_point-(adams_N-1)) to (end_point-1)
    */
void Integrator::GetDerivativeEnd(const std::vector<double>& f, std::vector<double>& df, int end_point)
{
    const double* R = lattice->R();

    // last two points just use quadratic formula
    for(int i = end_point - 2; i < end_point; i++)
    {
        double x = R[i];
        double a = R[i-2];
        double b = R[i-1];
        double c = R[i];

        double coeff_A = (x-b)/((a-b)*(a-c));
        double coeff_B = (x-a)/((b-a)*(b-c));
        double coeff_C = ((x-a) + (x-b))/((c-a)*(c-b));

        double A = f[i-2];
        double B = f[i-1];
        double C = f[i];

        df[i] = coeff_A * A + coeff_B * B + coeff_C * C;
        df[i] = df[i] * lattice->dR(i);
    }
}
