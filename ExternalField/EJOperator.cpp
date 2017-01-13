#include "EJOperator.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/MathConstant.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/factorials.hpp>

double EJOperator::GetReducedMatrixElement(const Orbital& b, const Orbital& a) const
{
    double matrix_element = 0.0;
    double coeff = 0.0;

    // This is the reduced matrix element of the spherical tensor CJ
    MathConstant* math = MathConstant::Instance();
    coeff = math->SphericalTensorReducedMatrixElement(b.Kappa(), a.Kappa(), K);

    // Convert from reduced matrix element to full matrix element
    //coeff *= pow(-1.0, e1.J() - e1.M()) * MathConstant::Instance()->Electron3j(e2.J(), e1.J(), J, e2.M(), -e1.M());

    // Don't bother computing the overlap if the angular part is zero
    if(coeff == 0)
    {
        return 0.0;
    }

    const double* R = integrator->GetLattice()->R();

    // p = omega/c, in atomic units
    double omega = fabs(a.Energy() - b.Energy());
    double p = omega/math->SpeedOfLightAU();

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
                bessel.f[i] = math->sph_bessel_small_limit(K, R[i]);
            }

            double overlap = integrator->GetPotentialMatrixElement(b, a, bessel);

            matrix_element = coeff * overlap * boost::math::double_factorial<double>(2. * K + 1.);
        }
        else
        {
            RadialFunction integrand(mmin(a.size(), b.size()));
            double kappa_diff = double(b.Kappa() - a.Kappa())/double(K + 1);

            for(unsigned int i = 0; i < integrand.size(); i++)
            {
                integrand.f[i] = (b.f[i] * a.f[i] + b.g[i] * a.g[i]) * boost::math::sph_bessel(K, R[i] * p)
                    + ((b.f[i] * a.g[i] + b.g[i] * a.f[i]) * kappa_diff + (b.f[i] * a.g[i] - b.g[i] * a.f[i])) * boost::math::sph_bessel(K + 1, R[i] * p);
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
            double kappa_diff = double(b.Kappa() - a.Kappa())/double(K + 1);

            for(unsigned int i = 0; i < integrand.size(); i++)
            {
                double pR = p * R[i];
                integrand.f[i] = - kappa_diff * (math->sph_bessel_prime(K, pR) + math->sph_bessel(K, pR)/pR) * (b.f[i] * a.g[i] + b.g[i] * a.f[i])
                    + (b.f[i] * a.g[i] - b.g[i] * a.f[i]) * K * math->sph_bessel(K, pR)/pR;
            }


            double overlap = integrator->Integrate(integrand);
            matrix_element = coeff * overlap * boost::math::double_factorial<double>(2 * K + 1)/pow(p, K);
        }
    }

    return matrix_element;
}

double MJOperator::GetReducedMatrixElement(const Orbital& b, const Orbital& a) const
{
    double matrix_element = 0.0;
    double coeff = 0.0;

    MathConstant* math = MathConstant::Instance();
    coeff = math->SphericalTensorReducedMatrixElement(-b.Kappa(), a.Kappa(), K);

    // Don't bother computing the overlap if the angular part is zero
    if(coeff == 0)
    {
        return 0.0;
    }

    const double* R = integrator->GetLattice()->R();

    // p = omega/c, in atomic units
    double omega = fabs(a.Energy() - b.Energy());
    double p = omega/math->SpeedOfLightAU();

    if(p < 1.e-6)
    {
        // For the same state to itself, a special case is needed to compute the
        // overlap, because dependence of the spherical bessel function on k
        // exactly cancels the k dependence added by reintroducing dimensionality
        RadialFunction integrand(mmin(b.size(), a.size()));
        double kappa_plus = double(b.Kappa() + a.Kappa())/double(K + 1);

        for(unsigned int i = 0; i < integrand.size(); i++)
        {
            integrand.f[i] = (b.f[i] * a.g[i] + b.g[i] * a.f[i]) * math->sph_bessel_small_limit(K, R[i]) * kappa_plus;
        }

        double overlap = integrator->Integrate(integrand);
        matrix_element = coeff * overlap * boost::math::double_factorial<double>(2 * K + 1);

        // This is the matrix operator with dimensions in atomic units
        matrix_element *= 2.0 * math->SpeedOfLightAU();
    }
    else
    {
        RadialFunction integrand(mmin(b.size(), a.size()));
        double kappa_plus = double(b.Kappa() + a.Kappa())/double(K + 1);

        for(unsigned int i = 0; i < integrand.size(); i++)
        {
            integrand.f[i] = (b.f[i] * a.g[i] + b.g[i] * a.f[i]) * boost::math::sph_bessel(K, p * R[i]) * kappa_plus;
        }

        double overlap = integrator->Integrate(integrand);
        matrix_element = coeff * overlap * boost::math::double_factorial<double>(2 * K + 1)/pow(p, K);

        // This is the matrix operator with dimensions in atomic units
        matrix_element *= 2.0 * math->SpeedOfLightAU();
    }

    return matrix_element;
}

EMCalculator::EMCalculator(MultipolarityType type, int J, MultirunOptions& user_input, pOrbitalManagerConst orbitals, pLevelStore levels, pIntegrator integrator):
    TransitionCalculator(user_input, orbitals, levels), type(type), J(J)
{
    if(type == MultipolarityType::E)
        op = std::make_shared<EJOperator>(J, integrator);
    else
        op = std::make_shared<MJOperator>(J, integrator);
}

void EMCalculator::PrintHeader() const
{
    // Print the appropriate header for the type of transitions we're printing:
    // either reduced matrix elements of transition line strengths (the
    // default)
    if(user_input.search("--reduced-elements")){
      *outstream << Name(type) << J << " reduced matrix elements (T):" 
          << std::endl;
    } else {
      *outstream << Name(type) << J << " transition strengths (S):" 
          << std::endl;
    }
}

void EMCalculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    MathConstant* math = MathConstant::Instance();

    int twoj1 = left.first->GetTwoJ();
    int twoj2 = right.first->GetTwoJ();
    double value = matrix_element/math->Electron3j(twoj2, twoj1, J, twoj2, -twoj1);

    // Read over the user input to see if we need to print the reduced matrix
    // elements or the transition line strengths
    if(user_input.search("--reduced-elements")){ 
      *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(6) << value << std::endl;
    } else {
      *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(6) << value * value << std::endl;

    }
}
