#include "Interpolator.h"
#include "Include.h"

Interpolator::Interpolator(const std::vector<double>& R_points, unsigned int order):
    lattice(NULL)
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
    R(R_points), dR(dR_points), lattice(NULL)
{}

void Interpolator::Interpolate(const std::vector<double>& yvalues, double xvalue, double& yvalue, double& derivative, unsigned int order)
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
    {   left = 0;
        if(xvalue > R[R.size()-1])
            left = R.size();
        else
        {   while(xvalue > R[left])
                left++;
        }
    }
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

    std::vector<double> f(order);
    std::vector<double> df(order);
    std::vector<double> x(order);

    unsigned int k;
    for(k = 0; k < order; k++)
    {   f[k] = yvalues[start_point+k];
        df[k] = 0.;
        if(lattice)
            x[k] = lattice->R(start_point+k);
        else
            x[k] = R[start_point+k];
    }

    std::vector<double> fnext(order);
    std::vector<double> dfnext(order);

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

void Interpolator::GetDerivative(const std::vector<double>& y, std::vector<double>& dydr, unsigned int order)
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
