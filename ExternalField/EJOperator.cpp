#include "EJOperator.h"
#include "Include.h"
#include "Universal/Function.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/MathConstant.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/factorials.hpp>

double EJOperator::GetMatrixElement(const Orbital& b, const SingleParticleWavefunction& a) const
{
    double matrix_element = 0.0;
    double coeff = 0.0;

    // This is the reduced matrix element of the spherical tensor CJ
    coeff = SphericalTensorReducedMatrixElement(K, b.Kappa(), a.Kappa());

    // Convert from reduced matrix element to full matrix element
    //coeff *= pow(-1.0, e1.J() - e1.M()) * MathConstant::Instance()->Electron3j(e2.J(), e1.J(), J, e2.M(), -e1.M());

    // Don't bother computing the overlap if the angular part is zero
    if(coeff == 0)
    {
        return 0.0;
    }

    const double alpha = constants->GetAlpha();
    const double* R = integrator->GetLattice()->R();

    // p = omega/c, in atomic units
    double omega = fabs(a.Energy() - b.Energy());
    double p = omega * alpha;

    if(gauge == TransitionGauge::Length)
    {
        if(p < 1.e-6)
        {
            // For the same state to itself, a special case is needed to compute the
            // overlap, because dependence of the spherical bessel function on k
            // exactly cancels the k dependence added by reintroducing dimensionality
            RadialFunction bessel(mmin(a.size(), b.size()));
            for(unsigned int i = 0; i < bessel.size(); i++)
            {
                bessel.f[i] = sph_bessel_small_limit(K, R[i]);
            }

            double overlap = integrator->GetPotentialMatrixElement(a, b, bessel);

            matrix_element = coeff * overlap * boost::math::double_factorial<double>(2. * K + 1.);
        }
        else
        {
            RadialFunction integrand(mmin(a.size(), b.size()));
            double kappa_diff = double(a.Kappa() - b.Kappa())/double(K + 1);

            for(unsigned int i = 0; i < integrand.size(); i++)
            {
                integrand.f[i] = (a.f[i] * b.f[i] + a.g[i] * b.g[i]) * boost::math::sph_bessel(K, R[i] * p)
                    + ((a.f[i] * b.g[i] + a.g[i] * b.f[i]) * kappa_diff + (a.f[i] * b.g[i] - a.g[i] * b.f[i])) * boost::math::sph_bessel(K + 1, R[i] * p);
            }

            double overlap = integrator->Integrate(integrand);
            matrix_element = coeff * overlap * boost::math::double_factorial<double>(2. * K + 1.)/pow(p, K);
        }
    }
    else if(gauge == TransitionGauge::Velocity)
    {
        if(p < 1.e-6)
        {
            // For the same state to itself, a special case is needed to compute the
            // overlap, because dependence of the spherical bessel function on k
            // exactly cancels the k dependence added by reintroducing dimensionality
            for(unsigned int x=0; x<mmin(a.size(), b.size()); x++)
            {
                //overlap += ((p1.f[x] * p2.f[x]) + (p1.g[x] * p2.g[x])) * sph_bessel_small_limit(J, R[x]) * dR[x];
            }

            matrix_element = ((double) boost::math::double_factorial<double>((2.0 * ((double) K)) + 1.0)) * matrix_element;
        }
        else
        {
            RadialFunction integrand(mmin(a.size(), b.size()));
            double kappa_diff = double(a.Kappa() - b.Kappa())/double(K + 1);

            for(unsigned int i = 0; i < integrand.size(); i++)
            {
                double pR = p * R[i];
                integrand.f[i] = - kappa_diff * (sph_bessel_prime(K, pR) + boost::math::sph_bessel(K, pR)/pR) * (a.f[i] * b.g[i] + a.g[i] * b.f[i])
                    + (a.f[i] * b.g[i] - a.g[i] * b.f[i]) * K * boost::math::sph_bessel(K, pR)/pR;
            }


            double overlap = integrator->Integrate(integrand);
            matrix_element = coeff * overlap * boost::math::double_factorial<double>(2 * K + 1)/pow(p, K);
        }
    }

    return matrix_element;
}
