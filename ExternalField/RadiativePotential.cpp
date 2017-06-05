#include "RadiativePotential.h"
#include "Include.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h>
#include <boost/math/special_functions/expint.hpp>

// accuracies for GSL
#define accuracy_rel 1.e-10
#define accuracy_abs 1.e-20

UehlingDecorator::UehlingDecorator(pHFOperator wrapped_hf, double nuclear_rms_radius, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy)
{
    GenerateStepUehling(nuclear_rms_radius);
}

UehlingDecorator::UehlingDecorator(pHFOperator wrapped_hf, const RadialFunction& density, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy)
{
    GenerateUehling(density);
}

void UehlingDecorator::GenerateStepUehling(double nuclear_rms_radius)
{
    directPotential.Clear();
    directPotential.resize(lattice->size());

    MathConstant* math = MathConstant::Instance();
    double alpha = physicalConstant->GetAlpha();
    const double* R = lattice->R();
    const double* dR = lattice->dR();
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

    size_t sizelimit = 1000;
    double prefactor = alpha * Z/math->Pi();

    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(sizelimit);
    gsl_set_error_handler_off();
    double cumulative_error = 0.;

    // Get interior
    int i = 0;
    while(R[i] < nuclear_radius)
    {
        parameters.r = R[i];
        double integral;
        double abserr;

        int ierr = gsl_integration_qagiu(&integrand, 1.0, accuracy_abs, accuracy_rel, sizelimit, workspace, &integral, &abserr);
        if(ierr)
            cumulative_error += fabs(prefactor/parameters.r * abserr * dR[i]);
        directPotential.f[i] = prefactor/parameters.r * integral;

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

        int ierr = gsl_integration_qagiu(&integrand, 1.0, accuracy_abs, accuracy_rel, sizelimit, workspace, &integral, &abserr);
        if(ierr)
            cumulative_error += fabs(prefactor/parameters.r * abserr * dR[i]);
        directPotential.f[i] = prefactor/parameters.r * integral;

        if(fabs(directPotential.f[i]) < 1.e-15)
            break;
        i++;
    }

    directPotential.resize(mmin(i+1, lattice->size()));

    // Print uncertainty
    double V_integrated = integrator->Integrate(directPotential);
    if(fabs(cumulative_error/V_integrated) > accuracy_rel)
        *logstream << std::setprecision(6) << "UehlingDecorator::GenerateStepUehling() error: "
                   << cumulative_error << " out of " << V_integrated << std::endl;

    // Get derivative
    differentiator->GetDerivative(directPotential.f, directPotential.dfdr);
}

void UehlingDecorator::GenerateUehling(const RadialFunction& density)
{
    directPotential.Clear();
    directPotential.resize(lattice->size());

    MathConstant* math = MathConstant::Instance();
    double alpha = physicalConstant->GetAlpha();
    const double* R = lattice->R();
    const double* dR = lattice->dR();

    // Functions for inner integral (over t) with r, rp fixed
    struct integration_parameters_type
    {   double r;
        double rp;
        double alpha;
    } parameters;

    parameters.alpha = alpha;

    // Interior integrand: params = {r, rp, alpha}
    gsl_function integrand;
    integrand.params = &parameters;
    integrand.function = [](double t, void* parameters){
        const integration_parameters_type* params = static_cast<const integration_parameters_type*>(parameters);
        double t2 = t * t;
        double t3 = t2 * t;
        double t5 = t3 * t2;
        double rdiff = 2.0 * t * fabs(params->rp - params->r)/params->alpha;
        double rsum = 2.0 * t * (params->rp + params->r)/params->alpha;

        double part1 = std::sqrt(t2 - 1.0) * (1.0/t3 + 0.5/t5);
        double part2 = std::exp(-rdiff) - std::exp(-rsum);

        return part1 * part2;
    };

    size_t sizelimit = 1000;
    double prefactor = gsl_pow_2(alpha)/(6. * math->Pi());

    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(sizelimit);
    gsl_set_error_handler_off();
    double cumulative_error = 0.;

    unsigned int i_r = 0;
    while(i_r < lattice->size())
    {
        parameters.r = R[i_r];
        const double& r = parameters.r;
        double cumulative_error_p = 0.;

        // Get outer integrand, phi(rp), which will be integrated with the usual Integrator.
        RadialFunction phi(density.size());
        for(int p = 0; p < phi.size(); p++)
        {
            parameters.rp = R[p];

            // Get phi
            double integral;
            double abserr;

            int ierr = gsl_integration_qagiu(&integrand, 1.0, accuracy_abs, accuracy_rel, sizelimit, workspace, &integral, &abserr);
            phi.f[p] = density.f[p]/parameters.rp * integral;

            if(ierr)
                cumulative_error_p += fabs(density.f[p]/parameters.rp * abserr * dR[p]);
        }

        directPotential.f[i_r] = prefactor/r * integrator->Integrate(phi);
        cumulative_error += fabs(prefactor/r * cumulative_error_p * dR[i_r]);

        if(fabs(directPotential.f[i_r]) < 1.e-15)
            break;
        i_r++;
    }

    directPotential.resize(mmin(i_r+1, lattice->size()));

    // Print uncertainty
    double V_integrated = integrator->Integrate(directPotential);
    if(fabs(cumulative_error/V_integrated) > accuracy_rel)
        *logstream << std::setprecision(6) << "UehlingDecorator::GenerateUehling(density) error: "
                   << cumulative_error << " out of " << V_integrated << std::endl;

    // Get derivative
    differentiator->GetDerivative(directPotential.f, directPotential.dfdr);
}

MagneticSelfEnergyDecorator::MagneticSelfEnergyDecorator(pHFOperator wrapped_hf, double nuclear_rms_radius, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy)
{
    GenerateStepMagnetic(nuclear_rms_radius);
}

MagneticSelfEnergyDecorator::MagneticSelfEnergyDecorator(pHFOperator wrapped_hf, const RadialFunction& density, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy)
{
    GenerateMagnetic(density);
}

void MagneticSelfEnergyDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    wrapped->GetODEFunction(latticepoint, fg, w);

    if(latticepoint < magnetic.size())
    {
        w[0] += magnetic.f[latticepoint] * fg.f[latticepoint];
        w[1] -= magnetic.f[latticepoint] * fg.g[latticepoint];
    }
}
void MagneticSelfEnergyDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(latticepoint < magnetic.size())
    {
        w_f[0] += magnetic.f[latticepoint];
        w_g[1] -= magnetic.f[latticepoint];
    }
}
void MagneticSelfEnergyDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(latticepoint < magnetic.size())
    {
        jacobian[0][0] += magnetic.f[latticepoint];
        jacobian[1][1] -= magnetic.f[latticepoint];

        dwdr[0] += magnetic.dfdr[latticepoint] * fg.f[latticepoint];
        dwdr[1] -= magnetic.dfdr[latticepoint] * fg.g[latticepoint];
    }
}

SpinorFunction MagneticSelfEnergyDecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = wrapped->ApplyTo(a);

    ta.resize(mmax(ta.size(), magnetic.size()));

    double alpha = physicalConstant->GetAlpha();

    for(int i = 0; i < magnetic.size(); i++)
    {
        ta.f[i] -= magnetic.f[i] * a.g[i] / alpha;
        ta.dfdr[i] -= (magnetic.dfdr[i] * a.g[i] + magnetic.f[i] * a.dgdr[i]) / alpha;
        ta.g[i] -= magnetic.f[i] * a.f[i] / alpha;
        ta.dgdr[i] -= (magnetic.dfdr[i] * a.f[i] + magnetic.f[i] * a.dfdr[i]) / alpha;
    }

    return ta;
}

void MagneticSelfEnergyDecorator::GenerateStepMagnetic(double nuclear_rms_radius)
{
    magnetic.Clear();
    magnetic.resize(lattice->size());

    MathConstant* math = MathConstant::Instance();
    const double alpha = physicalConstant->GetAlpha();
    const double* R = lattice->R();
    const double* dR = lattice->dR();
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
        double ck = 2.0 * t * params->rnuc/params->alpha;
        double y = 2.0 * t * params->r/params->alpha;

        double part1 = std::sqrt(t2 - 1.0) * gsl_pow_3(ck);
        double part2 = 0.5 * (1. + ck) * ((1.-y) * std::exp(y-ck) - (1.+y) * std::exp(-ck-y));
        part2 = part2/gsl_pow_2(y);
        double part3 = y/3.;

        return (part2 + part3)/part1;
    };

    size_t sizelimit = 1000;
    double prefactor = 3. * Z * alpha/math->Pi();

    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(sizelimit);
    gsl_set_error_handler_off();
    double cumulative_error = 0.;

    // Get interior
    int i = 0;
    while(R[i] < nuclear_radius)
    {
        parameters.r = R[i];
        double integral;
        double abserr;

        int ierr = gsl_integration_qagiu(&integrand, 1.0, accuracy_abs, accuracy_rel, sizelimit, workspace, &integral, &abserr);
        if(ierr)
            cumulative_error += fabs(prefactor * abserr * dR[i]);
        magnetic.f[i] = prefactor * integral;

        i++;
    }

    // Exterior integrand
    integrand.function = [](double t, void* parameters){
        const integration_parameters_type* params = static_cast<const integration_parameters_type*>(parameters);
        double t2 = t * t;
        double ck = 2.0 * t * params->rnuc/params->alpha;
        double y = 2.0 * t * params->r/params->alpha;

        double part1 = std::sqrt(t2 - 1.0) * gsl_pow_2(y);
        double part2 = 0.5 * (1. + y) * ((1.-ck) * std::exp(ck-y) - (1.+ck) * std::exp(-ck-y));
        part2 = part2/gsl_pow_3(ck);
        double part3 = 1./3.;

        return (part2 + part3)/part1;
    };

    // Get exterior
    while(i < lattice->size())
    {
        parameters.r = R[i];
        double integral;
        double abserr;

        int ierr = gsl_integration_qagiu(&integrand, 1.0, accuracy_abs, accuracy_rel, sizelimit, workspace, &integral, &abserr);
        if(ierr)
            cumulative_error += fabs(prefactor * abserr * dR[i]);
        magnetic.f[i] = prefactor * integral;

        if(fabs(magnetic.f[i]) < 1.e-15)
            break;
        i++;
    }

    magnetic.resize(mmin(i+1, lattice->size()));

    // Print uncertainty
    double V_integrated = integrator->Integrate(magnetic);
    if(fabs(cumulative_error/V_integrated) > accuracy_rel)
        *logstream << std::setprecision(6) << "MagneticSelfEnergyDecorator::GenerateStepMagnetic() error: "
                   << cumulative_error << " out of " << V_integrated << std::endl;

    // Get derivative
    differentiator->GetDerivative(magnetic.f, magnetic.dfdr);
}

void MagneticSelfEnergyDecorator::GenerateMagnetic(const RadialFunction& density)
{
    magnetic.Clear();
    magnetic.resize(lattice->size());

    MathConstant* math = MathConstant::Instance();
    const double alpha = physicalConstant->GetAlpha();
    const double* R = lattice->R();
    const double* dR = lattice->dR();

    // Functions for inner integral (over t) with r, rp fixed
    struct inner_parameters_type
    {   double r;
        double rp;
        double alpha;
    } inner_parameters;

    inner_parameters.alpha = alpha;

    // Inner integrands (over t) for I(r) and H(r) = dI/dr.
    // H(r) = H1(r) for rp < r
    //    and H2(r) for r < rp
    // params = {r, rp, alpha}
    gsl_function inner_integrand_I;
    inner_integrand_I.params = &inner_parameters;
    inner_integrand_I.function = [](double t, void* parameters){
        const inner_parameters_type* params = static_cast<const inner_parameters_type*>(parameters);
        double t2 = t * t;
        double t3 = t2 * t;
        double rsum = 2.0 * t * (params->rp + params->r)/params->alpha;
        double rdiff = 2.0 * t * fabs(params->rp - params->r)/params->alpha;

        double part1 = std::sqrt(t2 - 1.0) * t3;
        double part2 = std::exp(-rdiff) - std::exp(-rsum) - rsum + rdiff;

        return part2/part1;
    };

    gsl_function inner_integrand_H;
    inner_integrand_H.params = &inner_parameters;
    inner_integrand_H.function = [](double t, void* parameters){
        const inner_parameters_type* params = static_cast<const inner_parameters_type*>(parameters);
        double t2 = t * t;
        double rsum = 2.0 * t * (params->rp + params->r)/params->alpha;
        double rdiff = 2.0 * t * fabs(params->rp - params->r)/params->alpha;    // r > rp

        double part1 = std::sqrt(t2 - 1.0) * t2;
        double part2;
        if(params->r > params->rp)
            part2 = 2.0/params->alpha * (-std::exp(-rdiff) + std::exp(-rsum));
        else
            part2 = 2.0/params->alpha * (std::exp(-rdiff) + std::exp(-rsum) - 2.0);

        return part2/part1;
    };

    size_t sizelimit = 1000;
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(sizelimit);
    gsl_set_error_handler_off();
    double cumulative_error = 0.;

    double prefactor = gsl_pow_4(alpha)/(16.0 * math->Pi());

    unsigned int i_r = 0;
    while(i_r < lattice->size())
    {
        inner_parameters.r = R[i_r];
        const double& r = inner_parameters.r;
        double cumulative_error_H = 0.;
        double cumulative_error_I = 0.;

        // Get outer integrands. I call them phi(rp) and chi(rp) for H and I, respectively.
        // These will be integrated with the usual Integrator.
        RadialFunction phi(density.size());
        RadialFunction chi(density.size());
        for(int p = 0; p < phi.size(); p++)
        {
            inner_parameters.rp = R[p];

            // Get phi
            double integral;
            double abserr;

            int ierrH = gsl_integration_qagiu(&inner_integrand_H, 1.0, accuracy_abs, accuracy_rel, sizelimit, workspace, &integral, &abserr);
            if(ierrH)
                cumulative_error_H += fabs(density.f[p]/inner_parameters.rp * abserr * dR[p]);
            phi.f[p] = density.f[p]/inner_parameters.rp * integral;

            int ierrI = gsl_integration_qagiu(&inner_integrand_I, 1.0, accuracy_abs, accuracy_rel, sizelimit, workspace, &integral, &abserr);
            if(ierrI)
                cumulative_error_I += fabs(density.f[p]/inner_parameters.rp * abserr * dR[p]);
            chi.f[p] = density.f[p]/inner_parameters.rp * integral;
        }

        magnetic.f[i_r] = prefactor/r * (integrator->Integrate(phi) - integrator->Integrate(chi)/r);
        cumulative_error += fabs(prefactor/r) * (cumulative_error_H + cumulative_error_I/r) * dR[i_r];

        if(fabs(magnetic.f[i_r]) < 1.e-15)
            break;
        i_r++;
    }

    magnetic.resize(mmin(i_r+1, lattice->size()));

    // Print uncertainty
    double V_integrated = integrator->Integrate(magnetic);
    if(fabs(cumulative_error/V_integrated) > accuracy_rel)
        *logstream << std::setprecision(6) << "MagneticSelfEnergyDecorator::GenerateMagnetic() error: "
                   << cumulative_error << " out of " << V_integrated << std::endl;

    // Get derivative
    differentiator->GetDerivative(magnetic.f, magnetic.dfdr);
}

ElectricSelfEnergyDecorator::ElectricSelfEnergyDecorator(pHFOperator wrapped_hf, double nuclear_rms_radius, bool integrate_offmass_term, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy), integrate_offmass_term(integrate_offmass_term)
{
    Initialize();
    GenerateStepEhigh(nuclear_rms_radius);
    GenerateStepElow(nuclear_rms_radius);
}

ElectricSelfEnergyDecorator::ElectricSelfEnergyDecorator(pHFOperator wrapped_hf, const RadialFunction& density, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy), integrate_offmass_term(true)
{
    Initialize();
    GenerateEhigh(density);
    GenerateElow(density);
}

void ElectricSelfEnergyDecorator::GenerateStepEhigh(double nuclear_rms_radius)
{
    directPotential.Clear();
    directPotential.resize(lattice->size());

    MathConstant* math = MathConstant::Instance();
    const double alpha = physicalConstant->GetAlpha();
    const double* R = lattice->R();
    const double* dR = lattice->dR();
    double nuclear_radius = std::sqrt(5./3.) * nuclear_rms_radius/math->BohrRadiusInFermi();

    struct integration_parameters_type
    {   double r;
        double rnuc;
        double alpha;
        double Z;
    } parameters;

    parameters.rnuc = nuclear_radius;
    parameters.alpha = alpha;
    parameters.Z = Z;

    auto offmass = [this](double r) {
        return 1./(1. + Ra/r);
    };

    // Integrand: params = {r, r_n (atomic units), alpha}
    gsl_function integrand;
    integrand.params = &parameters;
    integrand.function = [](double t, void* parameters){
        const integration_parameters_type* params = static_cast<const integration_parameters_type*>(parameters);
        double t2 = t * t;
        double t3 = t2 * t;
        double ck = 2.0 * t * params->rnuc/params->alpha;
        double rdiff = 2.0 * t * fabs(params->rnuc - params->r)/params->alpha;
        double rsum = 2.0 * t * (params->rnuc + params->r)/params->alpha;

        double part1 = t3 * std::sqrt(t2 - 1.0);
        double part2 = (1. - 0.5/t2) * (std::log(t2-1) + 4.*log(1./params->Z/params->alpha + 0.5)) - 1.5 + 1./t2;
        double part3;
        if(params->r >= params->rnuc)
            part3 = (ck - 1.) * std::exp(-rdiff) + (ck + 1.) * std::exp(-rsum);
        else
            part3 = 2.*params->r/params->rnuc * ck - (1. + ck) * (std::exp(-rdiff) - std::exp(-rsum));

        return (part2 * part3)/part1;
    };

    size_t sizelimit = 1000;
    double prefactor = -AfitSP * 3. * Z * gsl_pow_4(alpha) / (16. * math->Pi() * gsl_pow_3(nuclear_radius));

    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(sizelimit);
    gsl_set_error_handler_off();
    double cumulative_error = 0.;

    // Get potential
    int i = 0;
    while(R[i] < lattice->size())
    {
        parameters.r = R[i];
        double integral;
        double abserr;

        int ierr = gsl_integration_qagiu(&integrand, 1.0, accuracy_abs, accuracy_rel, sizelimit, workspace, &integral, &abserr);
        if(ierr)
            cumulative_error += fabs(prefactor/parameters.r * abserr * dR[i]);

        directPotential.f[i] = prefactor/parameters.r * integral;

        if(parameters.r > nuclear_radius && fabs(directPotential.f[i]) < 1.e-15)
            break;
        i++;
    }

    directPotential.resize(mmin(i+1, lattice->size()));

    // Add offmass correction
    if(integrate_offmass_term)
    {
        RadialFunction correction(lattice->size());
        prefactor = AfitSP * 6. * Z * alpha * Ra / (4. * math->Pi() * gsl_pow_3(nuclear_radius));

        // Integrand over t
        struct offmass_parameters_type
        {   double r;
            double rp;
            double ra;
            double alpha;
            double Z;
        } offmass_parameters;

        offmass_parameters.alpha = alpha;
        offmass_parameters.Z = Z;
        offmass_parameters.ra = Ra;

        gsl_function offmass_integrand;
        offmass_integrand.params = &offmass_parameters;
        offmass_integrand.function = [](double t, void* parameters){
            const offmass_parameters_type* params = static_cast<const offmass_parameters_type*>(parameters);
            double t2 = t * t;
            double rdiff = 2.0 * t * fabs(params->rp - params->r)/params->alpha;
            double rsum = 2.0 * t * (params->rp + params->r)/params->alpha;
            double ra = 2.0 * t * params->ra/params->alpha;

            double part1 = std::exp(ra)/std::sqrt(t2 - 1.0);
            double part2 = (1. - 0.5/t2) * (std::log(t2-1) + 4.*log(1./params->Z/params->alpha + 0.5)) - 1.5 + 1./t2;
            double part3 = gsl_sf_expint_E1(rdiff + ra) - gsl_sf_expint_E1(rsum + ra);

            return part1 * part2 * part3;
        };

        int i = 0;
        while(i < lattice->size())
        {
            offmass_parameters.r = R[i];
            double cumulative_error_p = 0.;

            // Get outer integrand (over rp): call this phi
            RadialFunction phi(lattice->real_to_lattice(nuclear_radius)+1);
            for(int p = 0; p < phi.size(); p++)
            {
                offmass_parameters.rp = R[p];
                double integral;
                double abserr;

                int ierr = gsl_integration_qagiu(&offmass_integrand, 1.0, accuracy_abs, accuracy_rel, sizelimit, workspace, &integral, &abserr);
                if(ierr)
                    cumulative_error_p += fabs(offmass_parameters.rp * abserr * dR[p]);

                phi.f[p] = offmass_parameters.rp * integral;
            }

            correction.f[i] = prefactor/offmass_parameters.r * integrator->Integrate(phi);
            cumulative_error += fabs(prefactor/offmass_parameters.r * cumulative_error_p * dR[i]);

            if(offmass_parameters.r > nuclear_radius && correction.f[i] < 1.e-15)
                break;
            i++;
        }

        correction.resize(mmin(i+1, lattice->size()));
        directPotential += correction;
    }
    else
    {   for(i = 0; i < directPotential.size(); i++)
            directPotential.f[i] *= offmass(R[i]);
    }

    // Print uncertainty
    double V_integrated = integrator->Integrate(directPotential);
    if(fabs(cumulative_error/V_integrated) > accuracy_rel)
        *logstream << std::setprecision(6) << "ElectricSelfEnergyDecorator::GenerateStepEhigh() error: "
                   << cumulative_error << " out of " << V_integrated << std::endl;

    // Get derivative
    differentiator->GetDerivative(directPotential.f, directPotential.dfdr);
}

void ElectricSelfEnergyDecorator::GenerateEhigh(const RadialFunction& density)
{
    directPotential.Clear();
    directPotential.resize(lattice->size());

    MathConstant* math = MathConstant::Instance();
    const double alpha = physicalConstant->GetAlpha();
    const double* R = lattice->R();
    const double* dR = lattice->R();

    // Functions for inner integral (over t) with r, rp fixed
    struct inner_parameters_type
    {   double r;
        double rp;
        double ra;
        double alpha;
        double Z;
    } inner_parameters;

    inner_parameters.alpha = alpha;
    inner_parameters.ra = Ra;
    inner_parameters.Z = Z;

    // Inner integrand: params = {r, rp, ra, alpha, Z}
    gsl_function inner_integrand;
    inner_integrand.params = &inner_parameters;
    inner_integrand.function = [](double x, void* parameters){
        const inner_parameters_type* params = static_cast<const inner_parameters_type*>(parameters);
        double t = std::sqrt(1.+x*x)/x; // Transform 1
        double t2 = t * t;
        double rsum = 2.0 * t * (params->rp + params->r)/params->alpha;
        double rdiff = 2.0 * t * fabs(params->rp - params->r)/params->alpha;
        double ra = 2.0 * t * params->ra/params->alpha;

        //double common1 = std::sqrt(t2 - 1.0);   // No transform
        double common1 = x * std::sqrt(x*x + 1.0);  // Transform 1
        double common2 = (1. - 0.5/t2) * (std::log(t2-1) + 4.*log(1./params->Z/params->alpha + 0.5)) - 1.5 + 1./t2;

        double part1 = std::exp(-rsum) - std::exp(-rdiff);
        part1 = params->alpha/t * part1;

        double part2 = std::exp(ra) * (gsl_sf_expint_E1(rdiff + ra) - gsl_sf_expint_E1(rsum + ra));
        part2 = 2. * params->ra * part2;

        return common2/common1 * (part1 + part2);
    };

    size_t sizelimit = 1000;
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(sizelimit);
    gsl_set_error_handler_off();
    double cumulative_error = 0.;

    double prefactor = AfitSP * alpha / (4.0 * math->Pi());

    unsigned int i_r = 0;
    while(i_r < lattice->size())
    {
        inner_parameters.r = R[i_r];
        const double& r = inner_parameters.r;
        double cumulative_error_p = 0.;

        // Get outer integrand phi(rp), to be integrated with the usual Integrator.
        RadialFunction phi(density.size());
        for(int p = 0; p < phi.size(); p++)
        {
            inner_parameters.rp = R[p];

            // Get phi
            double integral;
            double abserr;

            int ierr = gsl_integration_qagiu(&inner_integrand, 0.0, accuracy_abs, accuracy_rel, sizelimit, workspace, &integral, &abserr);
            if(ierr)
                cumulative_error_p += fabs(density.f[p]/inner_parameters.rp * abserr * dR[p]);

            phi.f[p] = density.f[p]/inner_parameters.rp * integral;
        }

        directPotential.f[i_r] = prefactor * integrator->Integrate(phi)/r;
        cumulative_error += fabs(prefactor * cumulative_error_p/r * dR[i_r]);

        if(fabs(directPotential.f[i_r]) < 1.e-15)
            break;
        i_r++;
    }

    directPotential.resize(mmin(i_r+1, lattice->size()));

    // Print uncertainty
    double V_integrated = integrator->Integrate(directPotential);
    if(fabs(cumulative_error/V_integrated) > accuracy_rel)
        *logstream << std::setprecision(6) << "ElectricSelfEnergyDecorator::GenerateEhigh() error: "
                   << cumulative_error << " out of " << V_integrated << std::endl;

    // Get derivative
    differentiator->GetDerivative(directPotential.f, directPotential.dfdr);
}

/** Generate low-frequency electric part of self-energy with nuclear_rms_radius (fm). */
void ElectricSelfEnergyDecorator::GenerateStepElow(double nuclear_rms_radius)
{
    // Use point-like here
    RadialFunction Elow(lattice->size());

    const double alpha = physicalConstant->GetAlpha();
    const double* R = lattice->R();

    double prefactor = gsl_pow_3(Z * alpha) * Z;

    int i = 0;
    while(i < lattice->size())
    {
        Elow.f[i] = -prefactor * std::exp(-Z * R[i]);
        Elow.dfdr[i] = -Z * Elow.f[i];

        if(fabs(Elow.f[i]) < 1.e-15)
            break;
        i++;
    }

    Elow.resize(mmin(i+1, lattice->size()));

    directPotential += Elow * BfitSP;
    potDWave += Elow * BfitD;
}

/** Generate low-frequency electric part from given density. */
void ElectricSelfEnergyDecorator::GenerateElow(const RadialFunction& density)
{
    // Use point-like here
    RadialFunction Elow(lattice->size());

    double alpha = physicalConstant->GetAlpha();
    const double* R = lattice->R();

    double prefactor = gsl_pow_3(alpha) * Z / 2.;

    int i = 0;
    while(i < lattice->size())
    {
        const double& r = R[i];

        // Integrand is phi
        RadialFunction phi(density.size());
        for(int p = 0; p < phi.size(); p++)
        {
            const double& rp = R[p];
            double rsum = Z * (r + rp);
            double rdiff = Z * fabs(r - rp);

            phi.f[p] = density.f[p]/rp * ((rdiff + 1.) * std::exp(-rdiff) - (rsum + 1.) * std::exp(-rsum));
        }

        Elow.f[i] = -prefactor * integrator->Integrate(phi)/r;

        if(fabs(Elow.f[i]) < 1.e-15)
            break;
        i++;
    }

    Elow.resize(mmin(i+1, lattice->size()));

    // Get derivative
    differentiator->GetDerivative(Elow.f, Elow.dfdr);

    directPotential += Elow * BfitSP;
    potDWave += Elow * BfitD;
}

RadialFunction ElectricSelfEnergyDecorator::GetDirectPotential() const
{
    RadialFunction ret = wrapped->GetDirectPotential();
    ret += directPotential;

    return ret;
}

void ElectricSelfEnergyDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    wrapped->GetODEFunction(latticepoint, fg, w);

    if(latticepoint < directPotential.size())
    {
        const double alpha = physicalConstant->GetAlpha();
        switch(fg.L())
        {
            case 0:
            case 1:
                w[0] += alpha * directPotential.f[latticepoint] * fg.g[latticepoint];
                w[1] -= alpha * directPotential.f[latticepoint] * fg.f[latticepoint];
                break;
            case 2:
                w[0] += alpha * potDWave.f[latticepoint] * fg.g[latticepoint];
                w[1] -= alpha * potDWave.f[latticepoint] * fg.f[latticepoint];
            default:
                break;
        }
    }
}

void ElectricSelfEnergyDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(latticepoint < directPotential.size())
    {
        const double alpha = physicalConstant->GetAlpha();
        switch(fg.L())
        {
            case 0:
            case 1:
                w_g[0] += alpha * directPotential.f[latticepoint];
                w_f[1] -= alpha * directPotential.f[latticepoint];
                break;
            case 2:
                w_g[0] += alpha * potDWave.f[latticepoint];
                w_f[1] -= alpha * potDWave.f[latticepoint];
            default:
                break;
        }
    }
}

void ElectricSelfEnergyDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(latticepoint < directPotential.size())
    {
        const double alpha = physicalConstant->GetAlpha();
        switch(fg.L())
        {
            case 0:
            case 1:
                jacobian[0][1] += alpha * directPotential.f[latticepoint];
                jacobian[1][0] -= alpha * directPotential.f[latticepoint];
                dwdr[0] += alpha * directPotential.dfdr[latticepoint] * fg.g[latticepoint];
                dwdr[1] -= alpha * directPotential.dfdr[latticepoint] * fg.f[latticepoint];
                break;
            case 2:
                jacobian[0][1] += alpha * potDWave.f[latticepoint];
                jacobian[1][0] -= alpha * potDWave.f[latticepoint];
                dwdr[0] += alpha * potDWave.dfdr[latticepoint] * fg.g[latticepoint];
                dwdr[1] -= alpha * potDWave.dfdr[latticepoint] * fg.f[latticepoint];
            default:
                break;
        }
    }
}
void ElectricSelfEnergyDecorator::EstimateOrbitalNearOrigin(unsigned int numpoints, SpinorFunction& s) const
{
    RadialFunction V(wrapped->GetDirectPotential());

    switch(s.L())
    {
        case 0:
        case 1:
            V += directPotential;
            break;
        case 2:
            V += potDWave;
            break;
        default:
            break;
    }

    const int start_point = 0;
    const double alpha = physicalConstant->GetAlpha();

    double correction = 0.;
    if(s.size() >= numpoints)
        correction = s.f[start_point];
    else
        s.resize(numpoints);

    unsigned int i;
    for(i=start_point; i<start_point+numpoints; i++)
    {   if(s.Kappa() < 0)
    {   s.f[i] = pow(lattice->R(i), -s.Kappa());
        s.g[i] = alpha * s.f[i] * lattice->R(i) * V.f[i] / (2 * s.Kappa() - 1);
        s.dfdr[i] = - s.Kappa() * s.f[i] / lattice->R(i);
        s.dgdr[i] = ( - s.Kappa() + 1.) * s.g[i] / lattice->R(i);
    }
    else
    {   s.g[i] = alpha * pow(lattice->R(i), s.Kappa());
        s.f[i] = s.g[i] * lattice->R(i) * alpha * V.f[i] / (2 * s.Kappa() + 1);
        s.dgdr[i] = s.Kappa() * s.g[i] / lattice->R(i);
        s.dfdr[i] = (s.Kappa() + 1.) * s.f[i] / lattice->R(i);
    }
    }

    // Determine an appropriate scaling to make the norm close to unit.
    if(correction)
        correction = correction/s.f[start_point];
    else
        correction = Z * Z;

    for(i=start_point; i<start_point+numpoints; i++)
    {   s.f[i] = s.f[i] * correction;
        s.g[i] = s.g[i] * correction;
        s.dfdr[i] = s.dfdr[i] * correction;
        s.dgdr[i] = s.dgdr[i] * correction;
    }
}

SpinorFunction ElectricSelfEnergyDecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = wrapped->ApplyTo(a);
    switch(a.L())
    {
        case 0:
        case 1:
            ta -= a * directPotential;
            break;
        case 2:
            ta -= a * potDWave;
        default:
            break;
    }

    return ta;
}

void ElectricSelfEnergyDecorator::Initialize()
{
    const double alpha = physicalConstant->GetAlpha();
    Ra = 0.07 * Z*Z * gsl_pow_3(alpha);

    double Zp = (Z - 80.) * alpha;
    double Zp2 = Zp * Zp;
    double Zp3 = Zp2 * Zp;
    double Zp4 = Zp2 * Zp2;
    AfitSP = 1.071 - 1.976 * Zp2 - 2.128 * Zp3 + 0.169 * Zp4;

    double Za = Z * alpha;
    double Za2 = Za * Za;
    BfitSP = 0.074 + 0.35 * Za;
    BfitD  = 0.056 + 0.050 * Za + 0.195 * Za2;
}

QEDCalculator::QEDCalculator(MultirunOptions& user_input, Atom& atom):
    TransitionCalculator(user_input, atom)
{
    pHFOperatorConst hf = atom.GetHFOperator();
    pHFOperator zero = std::make_shared<HFOperatorBase>(*hf);
    pNucleusDecorator nucleus = atom.GetNucleusDecorator();

    double nuclear_rms_radius = user_input("NuclearRMSRadius", -1.0);
    if(nuclear_rms_radius < 0.0)
    {   if(nucleus)
            nuclear_rms_radius = nucleus->CalculateNuclearRMSRadius();
        else
            nuclear_rms_radius = 0.0;
    }

    pHFOperator qed = zero;

    if(user_input.search("--uehling"))
    {
        if(nucleus && user_input.search("--use-nuclear-density"))
            qed = std::make_shared<UehlingDecorator>(qed, nucleus->GetNuclearDensity());
        else
            qed = std::make_shared<UehlingDecorator>(qed, nuclear_rms_radius);
    }

    if(user_input.search("--self-energy"))
    {
        if(nucleus && user_input.search("--use-nuclear-density"))
        {
            if(!user_input.search("--no-magnetic"))
                qed = std::make_shared<MagneticSelfEnergyDecorator>(qed, nucleus->GetNuclearDensity());
            if(!user_input.search("--no-electric"))
                qed = std::make_shared<ElectricSelfEnergyDecorator>(qed, nucleus->GetNuclearDensity());
        }
        else
        {
            if(!user_input.search("--no-magnetic"))
                qed = std::make_shared<MagneticSelfEnergyDecorator>(qed, nuclear_rms_radius);
            if(!user_input.search("--no-electric"))
                qed = std::make_shared<ElectricSelfEnergyDecorator>(qed, nuclear_rms_radius);
        }
    }

    if(user_input.search("--rpa"))
        op = MakeRPA(qed, hf, atom.GetHartreeY());
    else
        op = qed;
}

void QEDCalculator::PrintHeader() const
{
    *outstream << "QED shift in 1/cm: " << std::endl;
}

void QEDCalculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    double value = matrix_element * MathConstant::Instance()->HartreeEnergyInInvCm();

    *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(12) << value << std::endl;
}
