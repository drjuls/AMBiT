#include "BasisGenerator.h"
#include "BSplineGrid.h"
#include "Include.h"
#include "Spline.h"
#include "Universal/SpinorFunction.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/Eigensolver.h"
#include "HartreeFock/StateIntegrator.h"
#include "Atom/MultirunOptions.h"

enum SplineType {NotreDame, Reno, Vanderbilt};
typedef Orbital BSpline;
typedef pOrbital pBSpline;

pStateManager BasisGenerator::GenerateBSplines(const std::vector<int>& max_pqn)
{
    excited->Clear();

    if(!max_pqn.size())
        return excited;

    bool debug = DebugOptions.OutputHFExcited();

    // Get spline type and parameters
    SplineType spline_type = Reno;
    std::string spline_type_string = user_input("Basis/BSpline/SplineType", "Reno");
    if(spline_type_string.compare("Reno") == 0 || spline_type_string.compare("DKB") == 0)
        spline_type = Reno;
    else if(spline_type_string.compare("Vanderbilt") == 0)
        spline_type = Vanderbilt;
    else if(spline_type_string.compare("NotreDame") == 0 || spline_type_string.compare("Johnson") == 0)
        spline_type = NotreDame;

    int n = user_input("Basis/BSpline/N", 40);
    int k = user_input("Basis/BSpline/K", 7);
    double rmax = user_input("Basis/BSpline/Rmax", lattice->R(lattice->Size()-1));

    if(spline_type == Reno)
    {   // Fix k. k should be at least lmax + 3
        int lmax = max_pqn.size() - 1;
        if(k < lmax + 3)
        {   k = lmax + 3;
            *errstream << "\nBSplineBasis: Warning: Order of Bsplines (k) too small for maximum l = " << lmax << ".";
            *errstream << "\n    k changed to " << k << "." << std::endl;
        }
    }

    double dr0;

    Eigensolver E;
    const double alpha = PhysicalConstant::Instance()->GetAlpha();
    const double alphasquared = PhysicalConstant::Instance()->GetAlphaSquared();
    pOPIntegrator integrator = hf->GetOPIntegrator();

    for(int l = 0; l < max_pqn.size(); l++)
    {
        if(!max_pqn[l])
            continue;

    for(int kappa = - (l+1); kappa <= l; kappa += 2*l + 1)
    {
        if(kappa == 0)
            break;

        if(spline_type == NotreDame)
        {
            // Set up grid and splines
            switch(l)
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

        // Total number of splines, including zeroed ones.
        unsigned int num_splines = n;
        if(spline_type == Reno)
            num_splines = n + abs(kappa) + 1;
        else if(spline_type == NotreDame)
            num_splines = n;
        else if(spline_type == Vanderbilt)
            num_splines = n + 1;

        // Create splines
        BSplineGrid grid(num_splines, k, dr0, rmax);
        int nderiv = 3;     // Number of derivatives of spline (including zeroth)
        double* fspline_buf = new double[nderiv * MaximumK];

        const double* knots = grid.GetSplineKnots();

        // Consistency check: make sure there are enough lattice points between each spline segment
        unsigned int prev_knot_point = 0;
        for(unsigned int i = 0; i < num_splines; i++)
        {
            unsigned int curr_knot_point = lattice->real_to_lattice(knots[i]);
            // skip zeros
            if((curr_knot_point - prev_knot_point) && (curr_knot_point - prev_knot_point) <= k)
            {   
                // Print whole list and escape!
                *errstream << "\nBSplineBasis: constructing Kappa = " << kappa << std::endl;
                *errstream << "Warning: too few points in each spline segment. Increase lattice size." << std::endl;
                unsigned int prev = 0;
                for(unsigned int i = 0; i < num_splines; i++)
                {   unsigned int curr = lattice->real_to_lattice(knots[i]);
                    if(curr > 0)
                        *errstream << i << " " << curr << "  " << curr - prev << "  " << knots[i] << std::endl;
                    prev = curr;
                }
                break;
            }
            else
                prev_knot_point = curr_knot_point;
        }

        unsigned int point = 0;
        unsigned int left = k-1;

        // Get values of splines on HF lattice
        unsigned int lattice_size = lattice->real_to_lattice(rmax);
        const double* R_lattice = lattice->R();

        std::vector<pBSpline> B_lattice(num_splines*2);
        std::vector<pBSpline>::iterator it = B_lattice.begin();
        while(it != B_lattice.end())
            *it++ = pBSpline(new BSpline(kappa, 0, 0.0, lattice_size));

        point = 0;
        left = k-1;

        while(point < lattice_size)
        {
            double x = R_lattice[point];

            while(knots[left+1] < x)
                left++;

            // 'left' is an array marker, so add 1 for fortran arrays
            int leftplus1 = (int)left + 1;
            bsplvd_(knots, &k, &x, &leftplus1, fspline_buf, &nderiv);

            // Transfer spline values from fspline_buf
            for(unsigned int s = 0; s < k; s++)
            {
                BSpline& Bupper = *B_lattice[s + left + 1 - k];      // B_i, i < n
                BSpline& Blower = *B_lattice[num_splines + s + left + 1 - k];  // B_i, i > n

                Bupper.f[point] = fspline_buf[s];
                Bupper.dfdr[point] = fspline_buf[s + MaximumK];

                Blower.g[point] = fspline_buf[s];
                Blower.dgdr[point] = fspline_buf[s + MaximumK];

                if(spline_type == Reno || spline_type == Vanderbilt)
                {
                    // Calculate small component
                    Bupper.g[point] = (Bupper.dfdr[point] + kappa/x * Bupper.f[point]) * alpha/2.;
                    Bupper.dgdr[point] = (fspline_buf[s + 2*MaximumK] + kappa/x * Bupper.dfdr[point] - kappa/(x*x) * Bupper.f[point]) * alpha/2.;

                    Blower.f[point] = (Blower.dgdr[point] - kappa/x * Blower.g[point]) * alpha/2.;
                    Blower.dfdr[point] = (fspline_buf[s + 2*MaximumK] - kappa/x * Blower.dgdr[point] + kappa/(x*x) * Blower.g[point]) * alpha/2.;
                }
            }

            point++;
        }

        delete[] fspline_buf;

        if(spline_type == Reno)
        {
            // Drop the splines at zero and Rmax to satisfy boundary conditions
            std::vector<pBSpline>::iterator it;
            it = B_lattice.begin();
            B_lattice.erase(it, it + abs(kappa));
            it = B_lattice.begin() + n;
            B_lattice.erase(it, it + abs(kappa) + 1);        
            B_lattice.erase(B_lattice.end() - 1);
        }
        else if(spline_type == Vanderbilt)
        {
            std::vector<pBSpline>::iterator it;
            it = B_lattice.begin();
            B_lattice.erase(it);
            it = B_lattice.begin() + n;
            B_lattice.erase(it);
        }

        // Wavefunctions are expanded in terms of B-splines as
        //     f = Sum_i (c_i. B_i)
        //     g = Sum_i (d_i. B_i)
        // To make integrals b(i,j) = <Bi|Bj> and A(i,j) = <Bi|H|Bj>

        unsigned int n2 = B_lattice.size();      // Size of the matrices A and b
        double* b = new double[n2*n2];
        double* A = new double[n2*n2];
        double* eigenvalues = new double[n2];

        unsigned int j;
        memset(b, 0, n2*n2*sizeof(double));
        memset(A, 0, n2*n2*sizeof(double));
        memset(eigenvalues, 0, n2*sizeof(double));


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
            BSpline& Bj = *B_lattice[j];
            SpinorFunction hf_applied_to_Bj = hf->ApplyTo(Bj);

            for(i=j; i<n2; i++)
            {
                BSpline& Bi = *B_lattice[i];
                A[i * n2 + j] = A[j * n2 + i] = integrator->GetInnerProduct(Bi, hf_applied_to_Bj);
                b[i * n2 + j] = b[j * n2 + i] = integrator->GetInnerProduct(Bi, Bj);
            }
        }

        // Boundary conditions to remove spurious states
        if(spline_type == NotreDame)
        {
            // Account for boundary conditions on splines.
            // At r = 0, this is effectively a delta function to push spurious states to high energy.
            if(kappa < 0)
                A[0] += 1./alpha;
            else
                A[0] += 2./alphasquared;
            A[0*n2 + n] += 0.5;
            A[n*n2 + 0] += -0.5;

            // At r = Rmax, these boundary conditions force f(r) = g(r) (effective mass->infinity).
            A[(n-1)*n2 + (n-1)] += 0.5/alpha;
            A[(n-1)*n2 + (n2-1)] += -0.5;
            A[(n2-1)*n2 + (n-1)] += 0.5;
            A[(n2-1)*n2 + (n2-1)] += -0.5*alpha;
        }
        else if(spline_type == Vanderbilt)
        {
            // At r = Rmax, these boundary conditions force f(r) = g(r) (effective mass->infinity).
            A[(n-1)*n2 + (n-1)] += 0.5/alpha;
            A[(n-1)*n2 + (n2-1)] += -0.5;
            A[(n2-1)*n2 + (n-1)] += 0.5;
            A[(n2-1)*n2 + (n2-1)] += -0.5*alpha;
        }

        // Solve A*p[i] = energy*B*p[i]
        if(E.SolveMatrixEquation(A, b, eigenvalues, n2))
        {
            if(debug)
                *outstream << "kappa = " << kappa << std::endl;

            bool reorth = user_input.search("Basis/BSpline/--reorthogonalise");

            // Interpolate from B-splines to HF lattice and store
            unsigned int pqn = l + 1;
            i = n2/2;
            pOrbitalConst s;

            while((pqn <= max_pqn[l]) && (i < n2))
            {
                // Check whether it is in the core
                s = open_core->GetState(OrbitalInfo(pqn, kappa));
                if(s)
                {   if(debug)
                    {   double diff = fabs((s->Energy() - eigenvalues[i])/s->Energy());
                        *outstream << "  " << s->Name() << " en: " << std::setprecision(8) << eigenvalues[i]
                                   << "  deltaE: " << diff << std::endl;
                    }

                    if(!closed_core->GetState(OrbitalInfo(pqn, kappa)))
                    {   pOrbital s_copy(new Orbital(s));
                        *s_copy = *s;
                        excited->AddState(s_copy);
                    }

                    pqn++;
                }
                else
                {   // Construct the orbital by summing splines with coefficients
                    pOrbital ds = pOrbital(new Orbital(kappa, pqn, eigenvalues[i]));
                    ds->ReSize(lattice_size);

                    double* coeff = &A[i*n2];

                    for(j=0; j<n2; j++)
                    {
                        (*ds) += (*B_lattice[j]) * (*coeff);
                        coeff++;
                    }

                    // Remove spurious states
                    if(fabs(ds->Norm(lattice) - 1.) > 1.e-2)
                    {   if(debug)
                            *outstream << "  SingleParticleWavefunction removed: energy = " << ds->Energy()
                                       << "  norm = " << ds->Norm(lattice) << std::endl;
                    }
                    else
                    {   if(reorth)
                            Orthogonalise(ds);
                        excited->AddState(ds);

                        if(debug)
                        {   *outstream << "  " << ds->Name() << " en: " << std::setprecision(8) << ds->Energy()
                                       << " norm: " << ds->Norm(lattice) - 1. << std::endl;
                        }

                        pqn++;
                    }
                }
                i++;
            }
        }
        else
        {   *errstream << "BSplineBasis: SolveMatrixEquation failed" << std::endl;
            exit(1);
        }
        
        delete[] eigenvalues;
        delete[] A;
        delete[] b;
        }   // kappa loop
    }

    return excited;
}
