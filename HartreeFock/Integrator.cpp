#include "Integrator.h"
#include "Include.h"

double Integrator::GetInnerProduct(const SpinorFunction& a, const SpinorFunction& b) const
{
    RadialFunction integrand = a.GetDensity(b);
    return Integrate(integrand);
}

double Integrator::GetInnerProduct(const RadialFunction& a, const RadialFunction& b) const
{
    RadialFunction integrand = a * b;
    return Integrate(integrand);
}

double Integrator::GetNorm(const SpinorFunction& a) const
{
    RadialFunction integrand = a.GetDensity();
    return Integrate(integrand);
}

double Integrator::GetPotentialMatrixElement(const SpinorFunction& a, const SpinorFunction& b, const RadialFunction& V) const
{
    RadialFunction integrand = a.GetDensity(b);
    integrand *= V;
    return Integrate(integrand);
}

double Integrator::GetPotentialMatrixElement(const SpinorFunction& a, const RadialFunction& V) const
{
    RadialFunction integrand = a.GetDensity();
    integrand *= V;
    return Integrate(integrand);
}

double SimpsonsIntegrator::Integrate(const RadialFunction& integrand) const
{
    return Integrate(integrand.size(), [&](int i){ return (integrand.f[i]); });
}

/** < a | b > = Integral (f_a * f_b + g_a * g_b) dr */
double SimpsonsIntegrator::GetInnerProduct(const SpinorFunction& a, const SpinorFunction& b) const
{
    int size = mmin(a.size(), b.size());
    return Integrate(size, [&](int i){ return (a.f[i] * b.f[i] + a.g[i] * b.g[i]); });
}

/** < a | b > = Integral (f_a * f_b) dr */
double SimpsonsIntegrator::GetInnerProduct(const RadialFunction& a, const RadialFunction& b) const
{
    int size = mmin(a.size(), b.size());
    return Integrate(size, [&](int i){ return (a.f[i] * b.f[i]); });
}

/** < a | a > */
double SimpsonsIntegrator::GetNorm(const SpinorFunction& a) const
{
    return Integrate(a.size(), [&](int i){ return (a.f[i] * a.f[i] + a.g[i] * a.g[i]); });
}

/** < a | V | b > = Integral (f_a * f_b + g_a * g_b) * V(r) dr */
double SimpsonsIntegrator::GetPotentialMatrixElement(const SpinorFunction& a, const SpinorFunction& b, const RadialFunction& V) const
{
    int size = mmin(a.size(), b.size());
    size = mmin(size, V.size());
    return Integrate(size, [&](int i){ return (a.f[i] * b.f[i] + a.g[i] * b.g[i]) * V.f[i]; });
}

/** < a | V | a > */
double SimpsonsIntegrator::GetPotentialMatrixElement(const SpinorFunction& a, const RadialFunction& V) const
{
    int size = mmin(a.size(), V.size());
    return Integrate(size, [&](int i){ return (a.f[i] * a.f[i] + a.g[i] * a.g[i]) * V.f[i]; });
}
