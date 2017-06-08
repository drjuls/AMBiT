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

    // Don't bother computing the overlap if the angular part is zero
    if(coeff == 0)
        return 0.0;

    const double* R = integrator->GetLattice()->R();

    // p = omega/c, in atomic units
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
                bessel.f[i] = gsl_pow_int(R[i], K);
            }

            double overlap = integrator->GetPotentialMatrixElement(b, a, bessel);
            matrix_element = coeff * overlap;
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
            matrix_element = coeff * overlap * boost::math::double_factorial<double>(2 * K + 1)/gsl_pow_int(p, K);
        }
    }
    else if(gauge == TransitionGauge::Velocity)
    {
        if(p > 1.e-6)
        {
            RadialFunction integrand(mmin(a.size(), b.size()));
            double kappa_diff = - double(b.Kappa() - a.Kappa())/double(K + 1);

            for(unsigned int i = 0; i < integrand.size(); i++)
            {
                double pR = p * R[i];
                integrand.f[i] = kappa_diff * (math->sph_bessel_prime(K, pR) + math->sph_bessel(K, pR)/pR) * (b.f[i] * a.g[i] + b.g[i] * a.f[i])
                    + (b.f[i] * a.g[i] - b.g[i] * a.f[i]) * K * math->sph_bessel(K, pR)/pR;
            }


            double overlap = integrator->Integrate(integrand);
            matrix_element = coeff * overlap * boost::math::double_factorial<double>(2 * K + 1)/gsl_pow_int(p, K);
        }
    }

    return matrix_element;
}

SpinorFunction EJOperator::ReducedApplyTo(const SpinorFunction& a, int kappa_b) const
{
    SpinorFunction ret(kappa_b);
    double coeff = 0.0;

    // This is the reduced matrix element of the spherical tensor CJ
    MathConstant* math = MathConstant::Instance();
    coeff = math->SphericalTensorReducedMatrixElement(kappa_b, a.Kappa(), K);

    // Don't bother computing the overlap if the angular part is zero
    if(coeff == 0)
        return ret;

    const double* R = integrator->GetLattice()->R();

    // p = omega/c, in atomic units
    double p = omega/math->SpeedOfLightAU();

    if(gauge == TransitionGauge::Length)
    {
        if(p < 1.e-6)
        {
            // For the same state to itself, a special case is needed to compute the
            // overlap, because dependence of the spherical bessel function on k
            // exactly cancels the k dependence added by reintroducing dimensionality
            RadialFunction bessel(a.size());
            for(unsigned int i = 0; i < a.size(); i++)
            {
                bessel.f[i] = gsl_pow_int(R[i], K);
                bessel.dfdr[i] = K * gsl_pow_int(R[i], K-1);
            }

            ret = a * bessel * coeff;
            ret.SetKappa(kappa_b);
        }
        else
        {
            double kappa_diff = double(kappa_b - a.Kappa())/double(K + 1);
            double kappa_diff_plus_one = kappa_diff + 1.;
            double kappa_diff_minus_one = kappa_diff - 1.;

            ret.resize(a.size());
            for(unsigned int i = 0; i < a.size(); i++)
            {
                double pR = p * R[i];
                double jK = math->sph_bessel(K, pR);
                double djKdr = p * math->sph_bessel_prime(K, pR);

                double jK1 = math->sph_bessel(K + 1, pR);
                double djK1dr = p * math->sph_bessel_prime(K + 1, pR);

                ret.f[i] = a.f[i] * jK + a.g[i] * jK1 * kappa_diff_plus_one;
                ret.dfdr[i] = a.dfdr[i] * jK + a.f[i] * djKdr
                              + (a.dgdr[i]  * jK1 + a.g[i] * djK1dr) * kappa_diff_plus_one;

                ret.g[i] = a.g[i] * jK + a.f[i] * jK1 * kappa_diff_minus_one;
                ret.dgdr[i] = a.dgdr[i] * jK + a.g[i] * djKdr
                              + (a.dfdr[i]  * jK1 + a.f[i] * djK1dr) * kappa_diff_minus_one;
            }

            ret *= coeff * boost::math::double_factorial<double>(2 * K + 1)/gsl_pow_int(p, K);
        }
    }
    else if(gauge == TransitionGauge::Velocity)
    {
        if(p > 1.e-6)
        {
            double kappa_diff = - double(kappa_b - a.Kappa())/double(K + 1);

            ret.resize(a.size());
            for(unsigned int i = 0; i < ret.size(); i++)
            {
                double pR = p * R[i];
                double jK = math->sph_bessel(K, pR);
                double jK1 = math->sph_bessel(K+1, pR);
                double jK2 = math->sph_bessel(K+2, pR);

                // Get primes and derivatives by recursion for speed
                double jKp = K/pR * jK - jK1;                                   //!< j_K'
                double jKpp = K*(K-1)/(pR*pR) * jK - (2*K+1)/pR * jK1 + jK2;    //!< j_K''

                double coeff_g = kappa_diff * jKp + (kappa_diff + K) * jK/pR;
                double coeff_f = kappa_diff * jKp + (kappa_diff - K) * jK/pR;

                ret.f[i] = coeff_g * a.g[i];
                ret.g[i] = coeff_f * a.f[i];

                ret.dfdr[i] = coeff_g * a.dgdr[i]
                    + p * (kappa_diff * jKpp + (kappa_diff + K) * (jKp/pR - jK/(pR*pR))) * a.g[i];
                ret.dgdr[i] = coeff_f * a.dfdr[i]
                    + p * (kappa_diff * jKpp + (kappa_diff - K) * (jKp/pR - jK/(pR*pR))) * a.f[i];
            }

            ret *= coeff * boost::math::double_factorial<double>(2 * K + 1)/pow(p, K);
        }
    }

    return ret;
}

SpinorFunction EJOperator::ConjugateReducedApplyTo(const SpinorFunction& a, int kappa_b) const
{
    SpinorFunction ret = ReducedApplyTo(a, kappa_b);
    if(gauge == TransitionGauge::Velocity)
        ret *= -1.;

    return ret;
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
            integrand.f[i] = (b.f[i] * a.g[i] + b.g[i] * a.f[i]) * gsl_pow_int(R[i], K);
        }

        double overlap = integrator->Integrate(integrand);
        matrix_element = coeff * overlap * kappa_plus;
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
    }

    // Return the matrix operator with dimensions in atomic units
    return matrix_element * 2.0 * math->SpeedOfLightAU();
}

SpinorFunction MJOperator::ReducedApplyTo(const SpinorFunction& a, int kappa_b) const
{
    SpinorFunction ret(kappa_b);
    double coeff = 0.0;

    MathConstant* math = MathConstant::Instance();
    coeff = math->SphericalTensorReducedMatrixElement(-kappa_b, a.Kappa(), K);

    // Don't bother computing the overlap if the angular part is zero
    if(coeff == 0)
    {
        return 0.0;
    }

    const double* R = integrator->GetLattice()->R();

    // p = omega/c, in atomic units
    double p = omega/math->SpeedOfLightAU();

    if(p < 1.e-6)
    {
        // For the same state to itself, a special case is needed to compute the
        // overlap, because dependence of the spherical bessel function on k
        // exactly cancels the k dependence added by reintroducing dimensionality
        double kappa_plus = double(kappa_b + a.Kappa())/double(K + 1);

        ret.resize(a.size());
        for(unsigned int i = 0; i < ret.size(); i++)
        {
            double RK = gsl_pow_int(R[i], K);
            double dRKdr = K * gsl_pow_int(R[i], K-1);

            ret.f[i] = a.g[i] * RK;
            ret.dfdr[i] = a.dgdr[i] * RK + a.g[i] * dRKdr;
            ret.g[i] = a.f[i] * RK;
            ret.dgdr[i] = a.dfdr[i] * RK + a.f[i] * dRKdr;
        }

        ret *= coeff * kappa_plus * 2.0 * math->SpeedOfLightAU();
    }
    else
    {
        double kappa_plus = double(kappa_b + a.Kappa())/double(K + 1);

        ret.resize(a.size());
        for(unsigned int i = 0; i < ret.size(); i++)
        {
            double pR = p * R[i];
            double jK = math->sph_bessel(K, pR);
            double djKdr = p * math->sph_bessel_prime(K, pR);

            ret.f[i] = a.g[i] * jK;
            ret.dfdr[i] = a.dgdr[i] * jK + a.g[i] * djKdr;
            ret.g[i] = a.f[i] * jK;
            ret.dgdr[i] = a.dfdr[i] * jK + a.f[i] * djKdr;
        }

        ret *= coeff * kappa_plus * boost::math::double_factorial<double>(2 * K + 1)/pow(p, K)
                * 2.0 * math->SpeedOfLightAU();
    }

    return ret;
}

EMCalculator::EMCalculator(MultipolarityType type, int J, MultirunOptions& user_input, Atom& atom):
    TransitionCalculator(user_input, atom), type(type), J(J)
{
    pHFOperatorConst hf = atom.GetHFOperator();

    if(type == MultipolarityType::E)
        op = std::make_shared<EJOperator>(J, hf->GetIntegrator());
    else
        op = std::make_shared<MJOperator>(J, hf->GetIntegrator());

    frequency_dependent_op = true;
    if(user_input.search("--rpa"))
        op = MakeRPA(std::static_pointer_cast<TimeDependentSpinorOperator>(op), hf, atom.GetHartreeY());
    else
    {
        double omega = user_input("Frequency", -1.0);
        if(omega >= 0.0)
        {   std::static_pointer_cast<TimeDependentSpinorOperator>(op)->SetFrequency(omega);
            frequency_dependent_op = false;
        }
        else
        {   frequency_dependent_op = true;
        }
    }
}

void EMCalculator::PrintHeader() const
{
    // Print the appropriate header for the type of transitions we're printing:
    // either reduced matrix elements of transition line strengths (the
    // default)
    if(user_input.search("--reduced-elements")){
      *outstream << Name() << " reduced matrix elements (T):"
          << std::endl;
    } else {
      *outstream << Name() << " transition strengths (S):"
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
    if(user_input.search("--reduced-elements"))
    {
        *outstream << "  " << ::Name(left) << " -> " << ::Name(right)
                   << " = " << std::setprecision(6) << value << std::endl;
    }
    else
    {
        *outstream << "  " << ::Name(left) << " -> " << ::Name(right)
                   << " = " << std::setprecision(6) << value * value << std::endl;

    }
}

std::string EMCalculator::Name() const
{
    return ::Name(type) + itoa(J);
}
