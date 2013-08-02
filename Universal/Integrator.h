#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/math/special_functions.hpp>
#include "Lattice.h"
#include "Function.h"
#include "CoupledFunction.h"
#include "MathConstant.h"

class Integrator
{
    /** General integrating class for linear ODEs of the form
          dy/dx = f(y,x) = f.coeff1(x) * y + f.coeff2(x)
     */
public:
    Integrator(Lattice* lat, unsigned int adams_order = 10): lattice(lat)
    {   SetAdamsOrder(adams_order);
    }
    virtual ~Integrator(void) {}
    
    void SetAdamsOrder(unsigned int adams_order);
    unsigned int AdamsOrder() { return adams_N; }

    /** Integrate - dy/dx = f(x)
        if (end_point > start_point)
            Integrate from start_point to end_point with y initialised
            from start_point-(adams_N-1) to start_point-1.
        if (end_point < start_point) (backwards integration)
            Integrate from start_point to end_point with y initialised
            from start_point+(adams_N-1) to start_point+1.
     */
    virtual void Integrate0(const std::vector<double>& f, std::vector<double>& y,
                       int start_point, int end_point);

    /** Integrate - dy/dx = f.coeff1(x) * y + f.coeff2(x)
        if (end_point > start_point)
            Integrate from start_point to end_point with y, dy initialised
            from start_point-(adams_N-1) to start_point-1.
        if (end_point < start_point) (backwards integration)
            Integrate from start_point to end_point with y, dy initialised
            from start_point+(adams_N-1) to start_point+1.
     */
    virtual void Integrate(const Function2& f, std::vector<double>& y, std::vector<double>& dy,
                       int start_point, int end_point);

    /** Integrate coupled linear ODE's
            dF/dr = A1.F + A2.G + A3
            dG/dr = A4.F + A5.G + A6
        where A1 .. A6 are functions of r only.
        if(end_point > start_point)
            Integrate from start_point to end_point with s.f, s.df, s.g, s.dg initialised
            from start_point-(adams_N-1) to start_point-1.
        if (end_point < start_point) (backwards integration)
            Integrate from start_point to end_point with s.f, s.df, s.g, s.dg initialised
            from start_point+(adams_N-1) to start_point+1.
    */
    virtual void Integrate2(const Function6& A, CoupledFunction& s, int start_point, int end_point);

    /** Calculate derivative df of any function f from start_point up to and NOT including end_point.
        This requires f to be known from (start_point-2) to (end_point+1).
     */
    virtual void GetDerivative(const std::vector<double>& f, std::vector<double>& df, int start_point, int end_point);

    /** Calculate the first two points of the derivative of f:
            df[start_point], df[start_point+1].
        Requires f to be known from (start_point) to (start_point+3)
     */
    virtual void GetDerivativeStart(const std::vector<double>& f, std::vector<double>& df, int start_point);

    /** Calculate the last two points of the derivative of f:
            df[end_point-2], df[end_point-1].
        Requires f to be known from (end_point-(adams_N-1)) to (end_point-1)
     */
    virtual void GetDerivativeEnd(const std::vector<double>& f, std::vector<double>& df, int end_point);

    virtual double BracketIntegral(const CoupledFunction& s1, const CoupledFunction& s2, double (*function)(double), int start_point, int end_point);
    virtual double BracketIntegral(const CoupledFunction& s1, const CoupledFunction& s2, boost::function < double (double r) > function, int start_point, int end_point, double f1factor = 1.0, double f2factor = 1.0, double g1factor = 1.0, double g2factor = 1.0, double crossfactor = 0.0);
    virtual double BracketIntegral(const CoupledFunction& s1, const CoupledFunction& s2, boost::function < double (double r) > *function, int start_point, int end_point, double f1factor = 1.0, double f2factor = 1.0, double g1factor = 1.0, double g2factor = 1.0, double crossfactor = 0.0);

    double AIntegralUpper(const CoupledFunction& s1, const CoupledFunction& s2, unsigned int J, double k);
    double AIntegralLower(const CoupledFunction& s1, const CoupledFunction& s2, unsigned int J, double k);
    double BIntegralUpper(const CoupledFunction& s1, const CoupledFunction& s2, unsigned int J, double k);
    double BIntegralLower(const CoupledFunction& s1, const CoupledFunction& s2, unsigned int J, double k);
protected:
    Lattice* lattice;
    unsigned int adams_N;

    const double* adams_coeff;  // Size = adams_N
};

#endif
