#include "Interpolator.h"
#include "Include.h"

Interpolator::Interpolator(const std::vector<double>& R_points, unsigned int order):
    lattice(pLattice())
{
    // Get lattice derivative dR by using the Interpolate function
    // Set the base grid to linear mesh x, and interpolate the R_points
    // in order to get their derivative, dR/dx = dR.
    std::vector<double> x(R_points.size());
    unsigned int i;
    for(i = 0; i < R_points.size(); i++)
        x[i] = double(i);

    R = x;
    dR.resize(R.size());
    double rval, deriv;
    for(i = 0; i < R_points.size(); i++)
    {   Interpolate(R_points, i, rval, deriv, order);
        dR[i] = deriv;
    }

    R = R_points;
}

Interpolator::Interpolator(const std::vector<double>& R_points, const std::vector<double>& dR_points):
    R(R_points), dR(dR_points), lattice(pLattice())
{}

void Interpolator::Interpolate(const std::vector<double>& yvalues, double xvalue, double& yvalue, double& derivative, unsigned int order) const
{
    // This routine uses Aitken method of interpolation.
    // (Lagrange interpolating polynomial)
    if(order%2 == 1)
        order++;

    unsigned int left;
    // get first lattice point greater than or equal to xvalue
    if(lattice)
        left = lattice->real_to_lattice(xvalue);
    else
        left = std::lower_bound(R.begin(), R.end(), xvalue) - R.begin();

    // left is the lattice point just less than xvalue
    if(left > 0)
        left--;

    unsigned int mid = order/2;
    unsigned int start_point;
    if(mid > left)
        start_point = 0;
    else if(int(mid) >= int(yvalues.size()) - int(left))
        start_point = yvalues.size()-order;
    else
        start_point = left - mid + 1;

    static std::vector<double> f, df, x;
    f.resize(order);
    df.resize(order);
    x.resize(order);

    unsigned int k;
    for(k = 0; k < order; k++)
    {   f[k] = yvalues[start_point+k];
        df[k] = 0.;
        if(lattice)
            x[k] = lattice->R(start_point+k);
        else
            x[k] = R[start_point+k];
    }

    static std::vector<double> fnext, dfnext;
    fnext.resize(order);
    dfnext.resize(order);

    for(unsigned int j = 1; j < order; j++)
    {    
        for(k = j; k < order; k++)
        {
            double dx = (xvalue - x[j-1])/(x[k] - x[j-1]);
            fnext[k] = (1.0 - dx)*f[j-1] + dx*f[k];
            dfnext[k] = (1.0 - dx)*df[j-1] + dx*df[k] + (f[k] - f[j-1])/(x[k] - x[j-1]);
        }
        for(k = j; k < order; k++)
        {   f[k] = fnext[k];
            df[k] = dfnext[k];
        }
    }
    
    yvalue = f[order - 1];
    derivative = df[order - 1];
}

void Interpolator::GetDerivative(const std::vector<double>& y, std::vector<double>& dydr, unsigned int order) const
{
    if(lattice)
    {   const double* R = lattice->R();

        for(unsigned int i=0; i<y.size(); i++)
        {   double yvalue, dyvalue;
            Interpolate(y, R[i], yvalue, dyvalue, order);
            dydr[i] = dyvalue;
        }
    }
    else
    {   for(unsigned int i=0; i<y.size(); i++)
        {   double yvalue, dyvalue;
            Interpolate(y, R[i], yvalue, dyvalue, order);
            dydr[i] = dyvalue;
        }
    }
}

double Interpolator::FindExtremum(const std::vector<double>& function, unsigned int zero_point, double& x) const
{
    // TODO: change order
    double f1 = function[zero_point-1];
    double f2 = function[zero_point];
    double f3 = function[zero_point+1];
    double f4 = function[zero_point+2];

    double A = (f4 - 3*f3 + 3*f2 - f1)/(6 * pow(lattice->H(), 3.));
    double B = (f3 - 2*f2 + f1)/(2 * lattice->H() * lattice->H());
    double C = (-f4 + 6*f3 - 3*f2 - 2*f1)/(6 * lattice->H());
    double D = f2;

    double x_offset = 0.;
    if(fabs(B/D) > 1.e-10)
    {   if(fabs(A/D) > 1.e-10)
            x_offset = -B/(3.*A)*(1.-sqrt(1. - 3.*A*C/(B*B)));
        else
            x_offset = -C/(2*B);
    }
    x = lattice->R(zero_point) + x_offset;

    return A*x_offset*x_offset*x_offset + B*x_offset*x_offset + C*x_offset + D;
}

