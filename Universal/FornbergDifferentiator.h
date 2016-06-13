#ifndef FORNBERG_DIFFERENTIATOR_H
#define FORNBERG_DIFFERENTIATOR_H

#include "Lattice.h"
#include <Eigen/Eigen>

class FornbergDifferentiator : public LatticeObserver
{
public:
    FornbergDifferentiator(pLattice lattice, int order = 5, bool include_second_derivative = false):
        LatticeObserver(lattice), order(order)
    {
        InitialiseCoefficients(include_second_derivative);
    }

    /** Get derivative of f.
        PRE: f is defined on the lattice, and has size not larger than it.
             dfdr.size() == f.size().
     */
    void GetDerivative(const std::vector<double>& f, std::vector<double>& dfdr) const;

    /** Get second derivative of f.
        Note that for this to work, include_second_derivative has to be set when initialising.
     */
    void GetSecondDerivative(const std::vector<double>& f, std::vector<double>& d2fdr2) const;

    virtual void Alert();

protected:
    void InitialiseCoefficients(bool include_second_derivative);

    /** Compute coefficients for finite difference approximation for the
        derivative of order k at xbar based on grid values at points in x.

        PRE: x[0]..x[n-1] are the grid points to be used, with n > k.
             The x values must be distinct.
             C is a matrix of size n by (k+1).
        POST: The m^th column vector of C (dimension 1 by n) contains the coefficients
             to approximate the m^th derivative of u at xvalue, u^{(m)}(xvalue).
             based on n values of u at x[0]..x[n-1].
             If U is a vector containing u(x) at these n points, then C[m]*U
             will give the approximation to u^{(m)}(xvalue).

        Based on the program "weights" in
            B. Fornberg, "Calculation of weights in finite difference formulas", SIAM Review 40 (1998), pp. 685-691.
        See also
            http://faculty.washington.edu/rjl/fdmbook/matlab/fdcoeffF.m
        for a matlab implementation.
     */
    void GetWeights(int k, double xvalue, const double* x, Eigen::MatrixXd& C) const;

protected:
    int order;

    // Derivative weights
    Eigen::MatrixXd first_derivative;
    Eigen::MatrixXd second_derivative;
};

#endif
