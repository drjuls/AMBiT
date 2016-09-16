#include "BasisGenerator.h"
#include "BSplineBasis.h"
#include "Include.h"
#include "Universal/SpinorFunction.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"
#include "Atom/MultirunOptions.h"
#include <Eigen/Eigen>
#include <gsl/gsl_bspline.h>

// This file contains B-spline routines from BasisGenerator as well as BSplineBasis
pOrbitalMap BasisGenerator::GenerateBSplines(const std::vector<int>& max_pqn)
{
    pOrbitalMap excited(new OrbitalMap(lattice));

    if(!max_pqn.size())
        return excited;

    bool debug = DebugOptions.OutputHFExcited();

    // Get spline type and parameters
    SplineType spline_type = SplineType::Reno;
    std::string spline_type_string = user_input("Basis/BSpline/SplineType", "Reno");
    if(spline_type_string.compare("Reno") == 0 || spline_type_string.compare("DKB") == 0)
        spline_type = SplineType::Reno;
    else if(spline_type_string.compare("Vanderbilt") == 0)
        spline_type = SplineType::Vanderbilt;
    else if(spline_type_string.compare("NotreDame") == 0 || spline_type_string.compare("Johnson") == 0)
        spline_type = SplineType::NotreDame;

    int n = user_input("Basis/BSpline/N", 40);
    int k = user_input("Basis/BSpline/K", 7);
    double rmax = user_input("Basis/BSpline/Rmax", lattice->MaxRealDistance());
    double dr0  = user_input("Basis/BSpline/R0", 0.0);

    if(rmax > lattice->MaxRealDistance())
        lattice->resize(rmax);

    if(spline_type == SplineType::Reno)
    {   // Fix k. k should be at least lmax + 3
        int lmax = max_pqn.size() - 1;
        if(k < lmax + 3)
        {   k = lmax + 3;
            *errstream << "\nBSplineBasis: Warning: Order of Bsplines (k) too small for maximum l = " << lmax << ".";
            *errstream << "\n    k changed to " << k << "." << std::endl;
        }
    }

    bool reorth = user_input.search("Basis/BSpline/--reorthogonalise");
    pIntegrator integrator = hf->GetIntegrator();

    // Make splines and store
    BSplineBasis spline_maker(lattice, n, k, rmax, dr0, spline_type);

    for(int l = 0; l < max_pqn.size(); l++)
    {
        if(!max_pqn[l])
            continue;

        for(int kappa = - (l+1); kappa <= l; kappa += 2*l + 1)
        {
            if(kappa == 0)
                break;

            if(debug)
                *logstream << "kappa = " << kappa << std::endl;

            pOrbitalMap basis = spline_maker.GenerateBSplines(hf, kappa, max_pqn[l]);

            for(auto it = basis->begin(); it != basis->end(); it++)
            {
                pOrbital ds = it->second;

                // Check whether it is in the core
                pOrbitalConst s = open_core->GetState(it->first);
                if(s)
                {   if(debug)
                    {   double diff = fabs((s->Energy() - ds->Energy())/s->Energy());
                        *logstream << "  " << s->Name() << " en: " << std::setprecision(8) << ds->Energy()
                                   << "  deltaE: " << diff << std::endl;
                    }

                    if(!closed_core->GetOccupancy(it->first))
                    {   pOrbital s_copy = s->Clone();
                        excited->AddState(s_copy);
                    }
                }
                else
                {   if(reorth)
                        Orthogonalise(ds);

                    if(debug)
                    {   *logstream << "  " << ds->Name() << " en: " << std::setprecision(8) << ds->Energy()
                                   << " norm: " << ds->Norm(integrator) - 1. << std::endl;
                    }

                    ds->ReNormalise(integrator);

                    excited->AddState(ds);
                }
            }
        }
    }

    return excited;
}

void BSplineBasis::SetParameters(int n, int k, double rmax, SplineType method)
{
    SetParameters(n, k, rmax, 0.0, method);
}

void BSplineBasis::SetParameters(int n, int k, double rmax, double dr0, SplineType method)
{
    N = n;
    K = k;
    Rmax = rmax;
    dR0 = dr0;
    spline_type = method;
}

std::vector<pBSpline> BSplineBasis::MakeSplines(int kappa, pPhysicalConstant constants)
{
    // Splines to return
    std::vector<pBSpline> splines;

    // Number of splines to make. These will be trimmed at the end to get back to N.
    int n = N;

    if(spline_type == SplineType::Reno)
        n = N + abs(kappa) + 1;
    else if(spline_type == SplineType::NotreDame)
        n = N;
    else if(spline_type == SplineType::Vanderbilt)
        n = N + 1;

    // Get breakpoints (knot points)
    gsl_vector* breakpoints = gsl_vector_alloc(n - K + 2);

    double dr0 = dR0;
    if(dr0 < 1.e-10)
    {
        if(spline_type == SplineType::NotreDame)
        {
            // dr0 depends on L
            int L = (kappa > 0? kappa: -kappa-1);
            switch(L)
            {   case 0:
                case 1:
                    dr0 = 0.0001;
                    break;
                case 2:
                case 3:
                    dr0 = 0.001;
                    break;
                default:
                    dr0 = 0.01;
                    break;
            }
        }
        else
            dr0 = 0.0001;
    }

    // Get parameter beta and h in
    //     breakpoints[i] = beta * (exp(h*(i-k+1)) - 1)
    // Subject to constraints
    //     breakpoints[0] = 0
    //     breakpoints[1] = dr0
    //     breakpoints[last] = rmax
    double a = Rmax/dr0;
    double ipow = n - K + 1;

    double xold;
    double x = 0.5;
    do
    {   xold = x;
        x = xold - (1. + a * xold - gsl_pow_int(1. + xold, ipow))
                   /(a - ipow * gsl_pow_int(1. + xold, ipow-1));
    }while(fabs(x - xold) > 1.e-15);

    double h = log(1. + x);
    double beta = dr0/x;

    gsl_vector_set(breakpoints, 0, 0.0);

    for(unsigned int i = 1; i < breakpoints->size; i++)
        gsl_vector_set(breakpoints, i, beta * (exp(h*i) - 1.));

    // Consistency check: make sure there are enough lattice points between each breakpoint
    unsigned int prev_knot_point = 0;
    for(unsigned int i = 1; i < breakpoints->size; i++)
    {
        unsigned int curr_knot_point = lattice->real_to_lattice(gsl_vector_get(breakpoints, i));

        if((curr_knot_point - prev_knot_point) <= K)
        {
            // Print whole list and escape!
            *errstream << "\nBSplineBasis: constructing Kappa = " << kappa << std::endl;
            *errstream << "Warning: too few points in each spline segment. Increase Lattice/NumPoints." << std::endl;
            unsigned int prev = 0;
            for(unsigned int i = 0; i < breakpoints->size; i++)
            {   unsigned int curr = lattice->real_to_lattice(gsl_vector_get(breakpoints, i));
                if(curr > 0)
                    *errstream << i << " " << curr << "  " << curr - prev << "  " << gsl_vector_get(breakpoints, i) << std::endl;
                prev = curr;
            }
            break;
        }
        else
            prev_knot_point = curr_knot_point;
    }

    gsl_bspline_workspace* gsl_workspace = gsl_bspline_alloc(K, breakpoints->size);
    gsl_bspline_knots(breakpoints, gsl_workspace);

    gsl_vector_free(breakpoints);

    // Allocate splines
    splines.resize(n * 2);

    unsigned int lattice_size = lattice->real_to_lattice(Rmax);
    const double* R = lattice->R();

    auto it = splines.begin();
    while(it != splines.end())
    {   *it = std::make_shared<BSpline>(kappa, 0, 0.0, lattice_size);
        ++it;
    }

    // Fill splines
    int nderiv = 2;     // Get zeroth, first and second derivatives of splines

    gsl_matrix* dB = gsl_matrix_alloc(n, nderiv+1);

    const double alpha = constants->GetAlpha();

    for(int i = 0; i < lattice_size; i++)
    {
        const double& x = R[i];
        gsl_bspline_deriv_eval(x, nderiv, dB, gsl_workspace);

        for(int s = 0; s < n; s++)
        {
            // Get non-zero spline values
            double spline_value = gsl_matrix_get(dB, s, 0);
            if(spline_value)
            {
                BSpline& Bupper = *splines[s];      // B_i, i < n
                BSpline& Blower = *splines[s + n];  // B_i, i >= n

                double deriv_value = gsl_matrix_get(dB, s, 1);
                Bupper.f[i] = spline_value;
                Bupper.dfdr[i] = deriv_value;

                Blower.g[i] = spline_value;
                Blower.dgdr[i] = deriv_value;

                // Apply dual kinetic balance
                if(spline_type == SplineType::Reno || spline_type == SplineType::Vanderbilt)
                {
                    double second_deriv = gsl_matrix_get(dB, s, 2);

                    // Calculate small component
                    Bupper.g[i] = (Bupper.dfdr[i] + kappa/x * Bupper.f[i]) * alpha/2.;
                    Bupper.dgdr[i] = (second_deriv + kappa/x * Bupper.dfdr[i] - kappa/(x*x) * Bupper.f[i]) * alpha/2.;

                    Blower.f[i] = (Blower.dgdr[i] - kappa/x * Blower.g[i]) * alpha/2.;
                    Blower.dfdr[i] = (second_deriv - kappa/x * Blower.dgdr[i] + kappa/(x*x) * Blower.g[i]) * alpha/2.;
                }
            }
        }
    }

    if(spline_type == SplineType::Reno)
    {
        // Drop the splines at zero and Rmax to satisfy boundary conditions
        std::vector<pBSpline>::iterator it;
        it = splines.begin();
        splines.erase(it, it + abs(kappa));

        it = splines.begin() + N;
        splines.erase(it, it + abs(kappa) + 1);
        splines.erase(splines.end() - 1);
    }
    else if(spline_type == SplineType::Vanderbilt)
    {
        std::vector<pBSpline>::iterator it;
        it = splines.begin();
        splines.erase(it);

        it = splines.begin() + N;
        splines.erase(it);
    }

    gsl_matrix_free(dB);
    gsl_bspline_free(gsl_workspace);

    if(splines.size() != 2 * N)
    {   *errstream << "BSplineBasis::MakeSplines: Number of splines (" << splines.size() << ") does not equal 2N = " << 2 * N << std::endl;
        exit(1);
    }

    return splines;
}

// Generate B-splines from core states up to max_pqn.
// If max_pqn <= 0, return all states including negative energy Dirac sea.
pOrbitalMap BSplineBasis::GenerateBSplines(pHFOperatorConst hf, int kappa, int max_pqn)
{
    pOrbitalMap excited(new OrbitalMap(lattice));
    bool debug = DebugOptions.OutputHFExcited();

    pPhysicalConstant physical_constant = hf->GetPhysicalConstant();
    const double alpha = physical_constant->GetAlpha();
    const double alphasquared = physical_constant->GetAlphaSquared();

    pIntegrator integrator = hf->GetIntegrator();

    // Create splines
    std::vector<pBSpline> splines = MakeSplines(kappa, physical_constant);

    // Wavefunctions are expanded in terms of B-splines as
    //     f = Sum_i (c_i. B_i)
    //     g = Sum_i (d_i. B_i)
    // To make integrals b(i,j) = <Bi|Bj> and A(i,j) = <Bi|H|Bj>

    unsigned int n2 = splines.size();      // size of the matrices A and b
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n2, n2);
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(n2, n2);

    unsigned int j;

    /* jcb: I removed the Gaussian integration that was previously being used to calculate <Bi|Bj>
            as well as the direct part of the potential. I think using the normal lattice is good enough.
            If someone wants to restore the Gaussian integration, use something like
                double BB = B_gauss[i*size + point] * B_gauss[j*size + point] * dR_grid[point];
                int_1 += BB;
                int_V += BB * potential_grid[point];
                int_KappaOnR += BB * kappa/R_grid[point];
    */
    unsigned int i = 0;
    for(j=0; j<n2; j++)
    {
        BSpline& Bj = *splines[j];
        SpinorFunction hf_applied_to_Bj = hf->ApplyTo(Bj);

        for(i=j; i<n2; i++)
        {
            BSpline& Bi = *splines[i];
            A(i, j) = A(j, i) = integrator->GetInnerProduct(Bi, hf_applied_to_Bj);
            b(i, j) = b(j, i) = integrator->GetInnerProduct(Bi, Bj);
        }
    }

    // Boundary conditions to remove spurious states
    if(spline_type == SplineType::NotreDame)
    {
        // Account for boundary conditions on splines.
        // At r = 0, this is effectively a delta function to push spurious states to high energy.
        if(kappa < 0)
            A(0, 0) += 1./alpha;
        else
            A(0, 0) += 2./alphasquared;
        A(0, N) += 0.5;
        A(N, 0) += -0.5;

        // At r = Rmax, these boundary conditions force f(r) = g(r) (effective mass->infinity).
        A((N-1), (N-1)) += 0.5/alpha;
        A((N-1), (n2-1)) += -0.5;
        A((n2-1), (N-1)) += 0.5;
        A((n2-1), (n2-1)) += -0.5*alpha;
    }
    else if(spline_type == SplineType::Vanderbilt)
    {
        // At r = Rmax, these boundary conditions force f(r) = g(r) (effective mass->infinity).
        A((N-1), (N-1)) += 0.5/alpha;
        A((N-1), (n2-1)) += -0.5;
        A((n2-1), (N-1)) += 0.5;
        A((n2-1), (n2-1)) += -0.5*alpha;
    }

    // Solve A*p[i] = energy*B*p[i]
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A,b);
    const Eigen::VectorXd& eigenvalues = es.eigenvalues();
    const Eigen::MatrixXd& eigenvectors = es.eigenvectors();

    // Interpolate from B-splines to HF lattice and store
    int L = (kappa > 0? kappa: -kappa-1);
    int pqn = L + 1;
    i = N;

    // Include Dirac sea if max_pqn <= 0.
    if(max_pqn <= 0)
    {   pqn = pqn - N;
        i = 0;
    }

    while((i < n2) && ((max_pqn <= 0) || (pqn <= max_pqn)))
    {
        // Construct the orbital by summing splines with coefficients
        pOrbital ds = std::make_shared<Orbital>(kappa, pqn, eigenvalues[i]);

        for(j=0; j<n2; j++)
        {
            (*ds) += (*splines[j]) * eigenvectors(j, i);
        }
        
        // Remove spurious states
        if(i >= N && fabs(ds->Norm(integrator) - 1.) > 1.e-2)
        {   if(debug)
                *logstream << "  Orbital removed: kappa = " << kappa << "  energy = " << ds->Energy()
                           << "  norm = " << ds->Norm(integrator) << std::endl;
        }
        else
        {   excited->AddState(ds);
            pqn++;
        }
        i++;
    }

    return excited;
}

pOrbitalMap BSplineBasis::GeneratePositiveBasis(pHFOperatorConst hf, int kappa)
{
    return GenerateBSplines(hf, kappa, abs(kappa) + N);
}

pOrbitalMap BSplineBasis::GenerateCompleteBasis(pHFOperatorConst hf, int kappa)
{
    return GenerateBSplines(hf, kappa, 0);
}
