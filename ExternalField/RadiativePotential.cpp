#include "RadiativePotential.h"
#include "Include.h"
#include "Universal/Interpolator.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

UehlingDecorator::UehlingDecorator(pHFOperator wrapped_hf, double nuclear_rms_radius, pOPIntegrator integration_strategy):
    LocalPotentialDecorator<UehlingDecorator>(wrapped_hf, integration_strategy)
{
    GenerateStepUehling(nuclear_rms_radius);
}

void UehlingDecorator::GenerateStepUehling(double nuclear_rms_radius)
{
    directPotential.Clear();
    directPotential.resize(lattice->size());

    MathConstant* math = MathConstant::Instance();
    double alpha = physicalConstant->GetAlpha();
    const double* R = lattice->R();
    double nuclear_radius = std::sqrt(5./3.) * nuclear_rms_radius/math->BohrRadiusInFermi();

    struct integration_parameters_type
    {   double r;
        double rnuc;
        double alpha;
    } parameters;

    parameters.rnuc = nuclear_radius;
    parameters.alpha = alpha;

    // Interior integrand: params = {r, r_n (atomic units), alpha}
    gsl_function integrand;
    integrand.params = &parameters;
    integrand.function = [](double t, void* parameters){
        const integration_parameters_type* params = static_cast<const integration_parameters_type*>(parameters);
        double t2 = t * t;
        double t4 = t2 * t2;
        double ck = 2.0 * t * params->rnuc/params->alpha;

        double part1 = std::sqrt(t2 - 1.0) * (1.0/t2 + 0.5/t4);
        double part2 = 2.0/gsl_pow_3(ck);
        // Expanding sinh improves convergence.
        //      part3 = params->r/params->rnuc * ck - std::exp(-ck) * (1. + ck) * std::sinh(2. * t * params->r/params->alpha);
        double part3 = params->r/params->rnuc * ck - 0.5 * (1. + ck)
                        * (std::exp(2. * t * (params->r - params->rnuc)/params->alpha) - std::exp(-2. * t * (params->r + params->rnuc)/params->alpha));

        return part1 * part2 * part3;
    };

    double accuracy = 1.e-10;
    size_t sizelimit = 1000;
    double prefactor = alpha * Z/math->Pi();

    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(sizelimit);
    gsl_set_error_handler_off();

    // Get interior
    int i = 0;
    while(R[i] < nuclear_radius)
    {
        parameters.r = R[i];
        double integral;
        double abserr;

        int ierr = gsl_integration_qagiu(&integrand, 1.0, 0.0, accuracy, sizelimit, workspace, &integral, &abserr);
        if(ierr)
            *errstream << "UehlingDecorator: GSL QAGIU error: r = " << parameters.r << ", integral = " << integral << ", abserr = " << abserr << std::endl;
        directPotential.f[i] = prefactor/parameters.r * integral;

        *logstream << directPotential.f[i] << std::endl;
        i++;
    }

    // Exterior integrand
    integrand.function = [](double t, void* parameters){
        const integration_parameters_type* params = static_cast<const integration_parameters_type*>(parameters);
        double t2 = t * t;
        double t4 = t2 * t2;
        double ck = 2.0 * t * params->rnuc/params->alpha;

        double part1 = std::sqrt(t2 - 1.0) * (1.0/t2 + 0.5/t4);
        double part2 = 2.0/gsl_pow_3(ck);
        // Expanding sinh and cosh improves convergence
        //      part3 = std::exp(-2. * t * params->r/params->alpha) * (ck * std::cosh(ck) - std::sinh(ck));
        double part3 = 0.5 * ((ck - 1.0) * std::exp(2. * t * (params->rnuc - params->r)/params->alpha)
                            + (ck + 1.0) * std::exp(-2. * t * (params->rnuc + params->r)/params->alpha));

        return part1 * part2 * part3;
    };

    // Get exterior
    while(i < lattice->size())
    {
        parameters.r = R[i];
        double integral;
        double abserr;

        int ierr = gsl_integration_qagiu(&integrand, 1.0, 0.0, accuracy, sizelimit, workspace, &integral, &abserr);
        if(ierr)
            *errstream << "UehlingDecorator: GSL QAGIU error: r = " << parameters.r << ", integral = " << integral << ", abserr = " << abserr << std::endl;
        directPotential.f[i] = prefactor/parameters.r * integral;

        *logstream << directPotential.f[i] << std::endl;
        if(fabs(directPotential.f[i]) < 1.e-15)
            break;
        i++;
    }

    directPotential.resize(mmin(i+1, lattice->size()));

    // Get derivative
    Interpolator interp(lattice);
    interp.GetDerivative(directPotential.f, directPotential.dfdr, 6);
}
