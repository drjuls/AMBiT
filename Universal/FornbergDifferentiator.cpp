#include "FornbergDifferentiator.h"
#include "Include.h"

void FornbergDifferentiator::GetDerivative(const std::vector<double>& f, std::vector<double>& dfdr) const
{
    int size = mmin(f.size(), first_derivative.rows());
    Eigen::Map<const Eigen::VectorXd> mapped_f(f.data(), size);
    Eigen::Map<Eigen::VectorXd> mapped_dfdr(dfdr.data(), size);

    for(int i = 0; i < size; i++)
    {
        int leftmost = mmax(0, i - order/2);
        if(leftmost < 0)
            leftmost = 0;
        else if(leftmost + order > size)
            leftmost = size - order;

        mapped_dfdr[i] = first_derivative.row(i) * mapped_f.segment(leftmost, order);
    }
}

void FornbergDifferentiator::GetSecondDerivative(const std::vector<double>& f, std::vector<double>& d2fdr2) const
{
    if(second_derivative.rows() == 0)
    {   d2fdr2.clear();
        return;
    }

    int size = mmin(f.size(), second_derivative.rows());
    Eigen::Map<const Eigen::VectorXd> mapped_f(f.data(), size);
    Eigen::Map<Eigen::VectorXd> mapped_d2fdr2(d2fdr2.data(), size);

    for(int i = 0; i < size; i++)
    {
        int leftmost = mmax(0, i - order/2);
        if(leftmost < 0)
            leftmost = 0;
        else if(leftmost + order > size)
            leftmost = size - order;

        mapped_d2fdr2[i] = second_derivative.row(i) * mapped_f.segment(leftmost, order);
    }
}

void FornbergDifferentiator::InitialiseCoefficients(bool include_second_derivative)
{
    int k = (include_second_derivative? 2: 1);
    Eigen::MatrixXd C(order, k+1);

    int size = lattice->size();
    const double* R = lattice->R();

    first_derivative = Eigen::MatrixXd::Zero(size, order);
    if(include_second_derivative)
        second_derivative = Eigen::MatrixXd::Zero(size, order);
    else
        second_derivative = Eigen::MatrixXd::Zero(0, 0);

    // Get weights for each point in lattice
    for(int i = 0; i < size; i++)
    {
        int leftmost = mmax(0, i - order/2);
        if(leftmost < 0)
            leftmost = 0;
        else if(leftmost + order > size)
            leftmost = size - order;

        GetWeights(k, R[i], &R[leftmost], C);

        first_derivative.row(i) = C.col(1).transpose();
        if(include_second_derivative)
            second_derivative.row(i) = C.col(2).transpose();
    }
}

void FornbergDifferentiator::Alert()
{
    // Only need to react if lattice size has increased
    if(lattice->size() > first_derivative.rows())
    {
        bool include_second_derivative = (second_derivative.rows() > 0);
        InitialiseCoefficients(include_second_derivative);
    }
}

void FornbergDifferentiator::GetWeights(int k, double xvalue, const double* x, Eigen::MatrixXd& C) const
{
    int n = C.rows();

    double c1 = 1.0;
    double c4 = x[0] - xvalue;

    C.setZero();
    C(0, 0) = 1.0;
    for(int i = 1; i < n; i++)
    {
        int mn = mmin(i, k);
        double c2 = 1.0;
        double c5 = c4;
        c4 = x[i] - xvalue;

        for(int j = 0; j < i; j++)
        {
            double c3 = x[i] - x[j];
            c2 *= c3;

            if(j == i-1)
            {
                for(int s = mn; s > 0; s--)
                {
                    C(i, s) = c1 * (s*C(i-1, s-1) - c5*C(i-1, s))/c2;
                }
                C(i, 0) = -c1 * c5 * C(i-1, 0)/c2;
            }

            for(int s = mn; s > 0; s--)
            {
                C(j, s) = (c4*C(j, s) - s*C(j, s-1))/c3;
            }

            C(j, 0) = c4 * C(j, 0)/c3;
        }

        c1 = c2;
    }
}
