#include <cmath>
#include <cstring>
#include <sstream>
#include <boost/math/special_functions/bessel.hpp>
#include "Function.h"
#include "Constant.h"

Polynomial::Polynomial(std::map<double,double> CoefficientsAndIndices)
{
    AiPairs = CoefficientsAndIndices;
}

double Polynomial::Evaluate(double r)
{
    double value = 0;
    std::map<double,double>::iterator it = AiPairs.begin();
    
    while(it != AiPairs.end())
    {
        value += it->second * pow(r, it->first);
        it++;
    }
    
    return value;
}

std::string Polynomial::Print(char variable)
{
    std::stringstream ss;
    ss.flush();
    
    std::map<double,double>::iterator it = AiPairs.begin();
    std::map<double,double>::iterator temp_it;
    while(it != AiPairs.end())
    {
        if(it->second != 1.000 && it->second != 0.000)
        {
            ss << it->second << variable << "^" << it->first;
        }
        else if(it->second == 1.000)
        {
            ss << variable << "^" << it->first;
        }
        
        temp_it = it;
        temp_it++;
        if(temp_it != AiPairs.end() && temp_it->second > 0 && it->second != 0)
        {
            ss << "+";    
        }
        it++;
    }
    
    return ss.str();
}

double triangular_condition(double a, double b, double c)
{
    double result = 0.0;
    
    if(a >= fabs(b - c) && a <= (b + c))
        if(b >= fabs(c - a) && b <= (c + a))
            if(c >= fabs(a - b) && c <= (a + b))
                result = 1.0;
            
    return result;
}

double triangular_condition_even(double a, double b, double c)
{
    double result = triangular_condition(a, b, c);
    double sum = a + b + c;
    if(((unsigned int) (10 * sum))%10 != 0)
    {
        result = 0.0;
    }
    if(((unsigned int) sum)%2 != 0)
    {
        result = 0.0;
    }
    
    return result;
}

double sph_bessel(unsigned int v, double x)
{
    return boost::math::sph_bessel(v, x);  
}

// j_n'(z) = d/dz j_n(z) = (n/z) j_n(z) - j_n+1(z)
double sph_bessel_prime(unsigned int v, double x)
{
    return (((((double) v)/x) * boost::math::sph_bessel(v, x)) - boost::math::sph_bessel(v + 1, x));   
}

double AIntegralFunction(unsigned int J, double k, double r)
{
    if(k == 0)
    {
        return sph_bessel_small_limit(J - 1, r);
    }
    
    return pow(r, J) * ((boost::math::sph_bessel(J - 1, k * r)) - ((double(J)/double(J+1)) * (boost::math::sph_bessel(J + 1, k * r)))) / pow(k, J - 1);
}

double BIntegralFunction(unsigned int J, double k, double r)
{
    if(k == 0)
    {
        return sph_bessel_small_limit(J, r);
    }
    return pow(r, J) * boost::math::sph_bessel(J, k * r) / pow(k, J - 1);    
}

double SphericalTensorReducedMatrixElement(unsigned int rank, int kappa1, int kappa2)
{
    if(kappa1 == 0.0 || kappa2 == 0.0)
    {
        return 0.0;
    }
    
    unsigned int l1 = 0.0;
    unsigned int l2 = 0.0;
    double j1 = 0.0;
    double j2 = 0.0;
    
    if(kappa1 < 0)
    {
        j1 = double(-kappa1) - 0.5;
        l1 = -kappa1 - 1;
    }
    else
    {
        j1 = double(kappa1) - 0.5;
        l1 = kappa1;
    }
    
    if(kappa2 < 0)
    {
        j2 = double(-kappa2) - 0.5;
        l2 = -kappa2 - 1;
    }
    else
    {
        j2 = double(kappa2) - 0.5;
        l2 = kappa2;
    }
    
    if((l1 + rank + l2)%2 != 0)
    {
        return 0.0;
    }
    
    return pow(-1.0, j1 + 0.5) * pow(((2.0 * j1) + 1.0) * ((2.0 * j2) + 1.0), 0.5) * Constant::Wigner3j(j1, j2, double(rank), -0.5, 0.5, 0.0);
}

double delta_function(double i, double j)
{
    if(i == j)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

double sph_bessel_small_limit(unsigned int v, double x)
{
    return (pow(x, v) * (1.0/boost::math::double_factorial<double>((2 * v) + 1.0)));
}
