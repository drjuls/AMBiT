#ifndef FUNCTION_H
#define FUNCTION_H

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

#endif