#include "BSplineGrid.h"
#include "Include.h"
#include "Spline.h"

BSplineGrid::BSplineGrid(unsigned int n, unsigned int k, double dr0, double rmax)
{
    this->n = n;
    this->k = k;

    // Create Spline knots and gaussian coords
    t = new double[n+k];
    CreateSplineKnots(dr0, rmax);

    // Gaussian coordinates and weights on the interval (0, 1)
    xgauss = new double[k];
    wgauss = new double[k];

    int kk = (int)k;
    gauss(&kk, xgauss, wgauss);

    // BSpline coordinate grid
    unsigned int size = (n-k+1)*k;
    r = (double*)malloc(sizeof(double) * size);
    dr = (double*)malloc(sizeof(double) * size);

    unsigned int point = 0;
    for(unsigned int left = k-1; left < n; left++)
    {
        for(unsigned int m = 0; m < k; m++)
        {
            r[point] = t[left] + (t[left+1] - t[left])*xgauss[m];
            dr[point] = (t[left+1] - t[left])*wgauss[m];

            point++;
        }
    }
    NumPoints = size;
}

BSplineGrid::~BSplineGrid(void)
{
    delete[] t;
    delete[] xgauss;
    delete[] wgauss;
}

/** Calculate the value that r[i] should be. */
//double BSplineGrid::lattice_to_real(unsigned int i) const;

/** Calculate the lattice spacing at a point. */
//double BSplineGrid::calculate_dr(double r_point) const;

/** Resizes the lattice such that NumPoints > min_size. */
//void BSplineGrid::ReSize(unsigned int min_size);

void BSplineGrid::CreateSplineKnots(double dr0, double rmax)
{
    // Get parameter beta and h in
    //     t[i] = beta * (exp(h*(i-k+1)) - 1)
    // Subject to constraints
    //     t[k-1] = 0
    //     t[k] = dr0
    //     t[n] = rmax
    double a = rmax/dr0;
    double ipow = double(n-k+1);
    
    double xold;
    double x = 0.5;
    do
    {   xold = x;
        x = xold - (1. + a * xold - pow(1. + xold, ipow))
                   /(a - ipow * pow(1. + xold, ipow-1.));
    }while(fabs(x - xold) > 1.e-15);

    h = log(1. + x);
    beta = dr0/x;

    unsigned int i;
    for(i = 0; i < k; i++)
        t[i] = 0.;

    for(i = k; i < n; i++)
        t[i] = beta * (exp(h*(i-k+1)) - 1.);

    for(i = n; i < n+k-1; i++)
        t[i] = rmax;
}
