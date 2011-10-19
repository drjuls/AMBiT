#ifndef FUNCTION_H
#define FUNCTION_H

#include <map>
#include <vector>

class Function2
{
public:
    Function2(void) {}
    virtual ~Function2(void) {}

    virtual double Coeff1(int point) const = 0;
    virtual double Coeff2(int point) const = 0;

    virtual std::vector<double> Coeff1(void) const = 0;
    virtual std::vector<double> Coeff2(void) const = 0;
};

class Function6
{
public:
    Function6(void) {}
    virtual ~Function6(void) {}

    virtual double Coeff1(int point) const = 0;
    virtual double Coeff2(int point) const = 0;
    virtual double Coeff3(int point) const = 0;
    virtual double Coeff4(int point) const = 0;
    virtual double Coeff5(int point) const = 0;
    virtual double Coeff6(int point) const = 0;
};

class Polynomial
{
public:
    Polynomial(std::map<double,double> IndicesAndCoefficients);
    
    double Evaluate(double r);
    std::string Print(char variable = 'r');
protected:
    std::map<double,double> AiPairs;
};

inline double EvaluatePolynomial(Polynomial aPoly, double r)
{
    return aPoly.Evaluate(r);
}

inline double SimpleFunction(double r)
{
    return r;
}

double sph_bessel_prime(unsigned int v, double x);
double sph_bessel_small_limit(unsigned int v, double x);
double delta_function(double i, double j);

#endif
