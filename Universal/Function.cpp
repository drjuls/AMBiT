#include <cmath>
#include <cstring>
#include <sstream>
#include <boost/math/special_functions/bessel.hpp>
#include "Function.h"

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

double sph_bessel_prime(unsigned int v, double x)
{
    return (((((double) v)/x) * boost::math::sph_bessel(v, x)) - boost::math::sph_bessel(v + 1, x));   
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
    return (pow(x / 2.0, v) * (1.0/boost::math::double_factorial<double>((2 * v) + 1.0)));
}