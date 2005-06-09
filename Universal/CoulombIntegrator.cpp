#include "Include.h"
#include "CoulombIntegrator.h"

double CoulombIntegrator::CoulombFunction::Coeff1(int point) const
{
    if(FirstPass)
        return -(K + 1.)/lattice.R(point);
    else
        return K/lattice.R(point);
}

double CoulombIntegrator::CoulombFunction::Coeff2(int point) const
{
    double ret = (*Density)[point]/lattice.R(point);
    if(FirstPass)
        return ret;
    else
        return -ret;
}

std::vector<double> CoulombIntegrator::CoulombFunction::Coeff1(void) const
{
    std::vector<double> ret(Density->size());
    // Check that lattice is big enough
    lattice.R(Density->size() - 1);
    const double* R = lattice.R();

    if(FirstPass)
    {   for(unsigned int i=0; i<ret.size(); i++)
            ret[i] = -(K + 1.)/R[i];
    }
    else
    {   for(unsigned int i=0; i<ret.size(); i++)
            ret[i] = K/R[i];
    }
    return ret;
}

std::vector<double> CoulombIntegrator::CoulombFunction::Coeff2(void) const
{
    std::vector<double> ret(Density->size());
    // Check that lattice is big enough
    lattice.R(Density->size() - 1);
    const double* R = lattice.R();

    if(FirstPass)
    {   for(unsigned int i=0; i<ret.size(); i++)
            ret[i] = (*Density)[i]/R[i];
    }
    else
    {   for(unsigned int i=0; i<ret.size(); i++)
            ret[i] = -(*Density)[i]/R[i];
    }
    return ret;
}

void CoulombIntegrator::CoulombIntegrate(std::vector<double>& density, std::vector<double>& potential, unsigned int k, double charge)
{
    // First part of Coulomb integration
    CoulombFunction function(lattice, density, k);
    std::vector<double> y(density.size());
    std::vector<double> dy(density.size());
    Integrate(function, y, dy, (adams_N-1), density.size());

    // Renormalise
    if(charge != 0.)
    {   double norm = y[density.size()-1] * lattice.R(density.size()-1);
        norm = charge/norm;
        for(unsigned int i=0; i<density.size(); i++)
        {   y[i] = y[i] * norm;
            dy[i] = dy[i] * norm;
            density[i] = density[i] * norm;
        }
    }

    // Second part of Coulomb integration
    potential.clear();
    potential.resize(density.size());

    if(k == 0.)
    {   /* For k == 0 we can use a special trick to improve accuracy,
             dI/dr = -I1/r
         */
        for(unsigned int i=density.size()-1; i>density.size()-adams_N; i--)
        {   potential[i] = y[i];
            dy[i] = -y[i]/lattice.R(i)*lattice.dR(i);
        }
        function.SetDensity(y);
        function.SetDirection(false);
        Integrate(function, potential, dy, density.size()-adams_N, -1);
    }
    else
    {   dy.clear();
        dy.resize(density.size());
        function.SetDensity(density);
        function.SetDirection(false);
        unsigned int i;
        for(i=density.size()-1; i>density.size()-adams_N; i--)
            density[i] = 0.;
        Integrate(function, potential, dy, density.size()-adams_N, -1);
        for(i=0; i<density.size(); i++)
        {   potential[i] = potential[i] + y[i];
        }
    }
}

void CoulombIntegrator::FastCoulombIntegrate(const std::vector<double>& density, std::vector<double>& potential, unsigned int k)
{
    potential.clear();
    potential.resize(density.size());

    lattice.R(density.size()-1);
    const double* dR = lattice.dR();

    if(k != 0)
    {
        const double* R_k = lattice.Rpower(k);
        const double* R_k_1 = lattice.Rpower(k+1);

        double sum = 0., sum1 = 0.;
        int i, j;
        double f0 = 0., f1 = 0., f2, f3;
        for(i=0, j=1; i<(int)density.size()-1; i+=2)
        {
            j=i+1;
            f2 = density[i] * dR[i] * R_k[i] / 3.;
            f3 = density[j] * dR[j] * R_k[j] / 3.;
            sum = sum + f0 + 4.*f1 + f2;
            sum1 = sum1 + f1 + 4.*f2 + f3;
            potential[i] = potential[i] + sum/R_k_1[i];
            potential[j] = potential[j] + sum1/R_k_1[j];
            f0 = f2;
            f1 = f3;
        }

        sum = 0.; sum1 = 0.;
        f0 = 0.; f1 = 0.;
        for(i = density.size()-1; i>0; i-=2)
        {
            j=i-1;
            f2 = density[i] * dR[i] / R_k_1[i] / 3.;
            f3 = density[j] * dR[j] / R_k_1[j] / 3.;
            sum = sum + f0 + 4.*f1 + f2;
            sum1 = sum1 + f1 + 4.*f2 + f3;
            potential[i] = potential[i] + sum * R_k[i];
            potential[j] = potential[j] + sum1 * R_k[j];
            f0 = f2;
            f1 = f3;
        }
    }
    else
    {
        const double* R = lattice.R();

        double sum = 0., sum1 = 0.;
        int i, j;
        double f0 = 0., f1 = 0., f2, f3;
        for(i=0, j=1; i<(int)density.size()-1; i+=2)
        {
            j=i+1;
            f2 = density[i] * dR[i] / 3.;
            f3 = density[j] * dR[j] / 3.;
            sum = sum + f0 + 4.*f1 + f2;
            sum1 = sum1 + f1 + 4.*f2 + f3;
            potential[i] = potential[i] + sum/R[i];
            potential[j] = potential[j] + sum1/R[j];
            f0 = f2;
            f1 = f3;
        }

        sum = 0.; sum1 = 0.;
        f0 = 0.; f1 = 0.;
        for(i = density.size()-1; i>0; i-=2)
        {
            j=i-1;
            f2 = density[i] * dR[i] / R[i] / 3.;
            f3 = density[j] * dR[j] / R[j] / 3.;
            sum = sum + f0 + 4.*f1 + f2;
            sum1 = sum1 + f1 + 4.*f2 + f3;
            potential[i] = potential[i] + sum;
            potential[j] = potential[j] + sum1;
            f0 = f2;
            f1 = f3;
        }
    }
}

void CoulombIntegrator::FastCoulombIntegrate(const std::vector<double>& density, std::vector<double>& potential, unsigned int k, unsigned int density_size)
{
    const double* dR = lattice.dR();

    if(k != 0)
    {
        const double* R_k = lattice.Rpower(k);
        const double* R_k_1 = lattice.Rpower(k+1);

        double sum = 0., sum1 = 0.;
        int i, j;
        double f0 = 0., f1 = 0., f2, f3;
        i = 0; j = 1;
        while(i<(int)density_size-1)
        {
            f2 = density[i] * dR[i] * R_k[i] / 3.;
            f3 = density[j] * dR[j] * R_k[j] / 3.;
            sum = sum + f0 + 4.*f1 + f2;
            sum1 = sum1 + f1 + 4.*f2 + f3;
            potential[i] = sum/R_k_1[i];
            potential[j] = sum1/R_k_1[j];
            f0 = f2;
            f1 = f3;

            i+=2;
            j+=2;
        }
        // Tail
        sum = (sum + sum1)/2.;
        i = density_size;
        while(i<(int)potential.size())
        {   potential[i] = sum/R_k_1[i];
            i++;
        }

        sum = 0.; sum1 = 0.;
        f0 = 0.; f1 = 0.;
        i = (int)density_size-1;
        j = i-1;
        while(i>0)
        {
            f2 = density[i] * dR[i] / R_k_1[i] / 3.;
            f3 = density[j] * dR[j] / R_k_1[j] / 3.;
            sum = sum + f0 + 4.*f1 + f2;
            sum1 = sum1 + f1 + 4.*f2 + f3;
            potential[i] += sum * R_k[i];
            potential[j] += sum1 * R_k[j];
            f0 = f2;
            f1 = f3;

            i-=2;
            j-=2;
        }
    }
    else
    {
        const double* R = lattice.R();

        double sum = 0., sum1 = 0.;
        int i, j;
        double f0 = 0., f1 = 0., f2, f3;
        i = 0; j = 1;
        while(i<(int)density_size-1)
        {
            f2 = density[i] * dR[i] / 3.;
            f3 = density[j] * dR[j] / 3.;
            sum = sum + f0 + 4.*f1 + f2;
            sum1 = sum1 + f1 + 4.*f2 + f3;
            potential[i] = sum/R[i];
            potential[j] = sum1/R[j];
            f0 = f2;
            f1 = f3;

            i+=2;
            j+=2;
        }
        // Tail
        sum = (sum + sum1)/2.;
        i = (int)density_size;
        while(i<(int)potential.size())
        {   potential[i] = sum/R[i];
            i++;
        }

        sum = 0.; sum1 = 0.;
        f0 = 0.; f1 = 0.;
        i = (int)density_size-1;
        j = i - 1;
        while(i>0)
        {
            f2 = density[i] * dR[i] / R[i] / 3.;
            f3 = density[j] * dR[j] / R[j] / 3.;
            sum = sum + f0 + 4.*f1 + f2;
            sum1 = sum1 + f1 + 4.*f2 + f3;
            potential[i] += sum;
            potential[j] += sum1;
            f0 = f2;
            f1 = f3;

            i-=2;
            j-=2;
        }
    }
}
