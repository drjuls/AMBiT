#include "BSplineBasis.h"
#include "BSplineGrid.h"
#include "Include.h"
#include "Spline.h"
#include "Universal/Interpolator.h"
#include "Universal/CoupledFunction.h"
#include "Universal/Constant.h"
#include "Universal/Eigensolver.h"

void BSplineBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    if(!num_states_per_l.size())
        return;

    NumStatesPerL = num_states_per_l;
    bool debug = DebugOptions.OutputHFExcited();

    if(spline_type == Reno)
    {   // Fix k. k should be at least lmax + 3
        unsigned int lmax = num_states_per_l.size() - 1;
        if(k < lmax + 3)
        {   k = lmax + 3;
            *errstream << "\nBSplineBasis: Warning: Order of Bsplines (k) too small for maximum l = " << lmax << ".";
            *errstream << "\n    k changed to " << k << "." << std::endl;
        }
    }

    double dr0;

    Eigensolver E;

    for(unsigned int l=0; l<num_states_per_l.size(); l++)
    {
        if(!num_states_per_l[l])
            continue;

    for(int kappa = - int(l+1); kappa <= int(l); kappa += 2*int(l) + 1)
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
        unsigned int num_splines;
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
        const double* dR_lattice = lattice->dR();

        std::vector<BSpline> B_lattice(num_splines*2, BSpline(kappa, lattice_size));

        point = 0;
        left = k-1;

        std::vector<double> Potential(core->GetHFPotential());

        while(point < lattice_size)
        {
            double x = R_lattice[point];
            double dx = dR_lattice[point];

            while(knots[left+1] < x)
                left++;

            // 'left' is an array marker, so add 1 for fortran arrays
            int leftplus1 = (int)left + 1;
            bsplvd_(knots, &k, &x, &leftplus1, fspline_buf, &nderiv);

            // Transfer spline values from fspline_buf
            for(unsigned int s = 0; s < k; s++)
            {
                BSpline& Bupper = B_lattice[s + left + 1 - k];      // B_i, i < n
                BSpline& Blower = B_lattice[num_splines + s + left + 1 - k];  // B_i, i > n

                Bupper.f[point] = fspline_buf[s];
                Bupper.df[point] = fspline_buf[s + MaximumK];

                Blower.g[point] = fspline_buf[s];
                Blower.dg[point] = fspline_buf[s + MaximumK];

                if(spline_type == Reno || spline_type == Vanderbilt)
                {
                    // Calculate small component
                    Bupper.g[point] = (Bupper.df[point] + kappa/x * Bupper.f[point])/2.;
                    Bupper.dg[point] = (fspline_buf[s + 2*MaximumK] + kappa/x * Bupper.df[point] - kappa/(x*x) * Bupper.f[point])/2.;

                    Blower.f[point] = Constant::AlphaSquared * (Blower.dg[point] - kappa/x * Blower.g[point])/2.;
                    Blower.df[point] = Constant::AlphaSquared * (fspline_buf[s + 2*MaximumK] - kappa/x * Blower.dg[point] + kappa/(x*x) * Blower.g[point])/2.;
                }

                Blower.df[point] *= dx;
                Blower.dg[point] *= dx;
                Bupper.df[point] *= dx;
                Bupper.dg[point] *= dx;                
            }

            point++;
        }

        delete[] fspline_buf;

        if(spline_type == Reno)
        {
            // Drop the splines at zero and Rmax to satisfy boundary conditions
            std::vector<BSpline>::iterator it;
            it = B_lattice.begin();
            B_lattice.erase(it, it + abs(kappa));
            it = B_lattice.begin() + n;
            B_lattice.erase(it, it + abs(kappa) + 1);        
            B_lattice.erase(B_lattice.end() - 1);
        }
        else if(spline_type == Vanderbilt)
        {
            std::vector<BSpline>::iterator it;
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
            BSpline& Bj = B_lattice[j];

            CoupledFunction exchange;
            core->CalculateExchange(Bj, exchange);
            exchange.ReSize(Bj.Size());

            for(i=j; i<n2; i++)
            {
                BSpline& Bi = B_lattice[i];

                // Same as
                //      En = StateIntegrator::HamiltonianMatrixElement(Bi, Bj, *core);
                // but faster since the exchange term is only calculated once.

                double En = 0.;
                for(unsigned int p=0; p<mmin(Bi.Size(), Bj.Size()); p++)
                {
                    double EQ1 = (Bj.df[p]/dR_lattice[p] + kappa*Bj.f[p]/R_lattice[p])
                                - (2. + Constant::AlphaSquared*Potential[p])*Bj.g[p]
                                - Constant::AlphaSquared * exchange.g[p];
                    double EQ2 = (Bj.dg[p]/dR_lattice[p] - kappa*Bj.g[p]/R_lattice[p]) + Potential[p]*Bj.f[p] + exchange.f[p];

                    En = En - (Bi.f[p] * EQ2 - Bi.g[p] * EQ1)*dR_lattice[p];
                }

                A[i * n2 + j] = A[j * n2 + i] = En;
                b[i * n2 + j] = b[j * n2 + i] = Bi.Overlap(Bj, lattice);
            }
        }

        // Boundary conditions to remove spurious states
        if(spline_type == NotreDame)
        {
            // Account for boundary conditions on splines.
            // At r = 0, this is effectively a delta function to push spurious states to high energy.
            if(kappa < 0)
                A[0] += 1./Constant::Alpha;
            else
                A[0] += 2./Constant::AlphaSquared;
            A[0*n2 + n] += 0.5;
            A[n*n2 + 0] += -0.5;

            // At r = Rmax, these boundary conditions force f(r) = g(r) (effective mass->infinity).
            A[(n-1)*n2 + (n-1)] += 0.5/Constant::Alpha;
            A[(n-1)*n2 + (n2-1)] += -0.5;
            A[(n2-1)*n2 + (n-1)] += 0.5;
            A[(n2-1)*n2 + (n2-1)] += -0.5*Constant::Alpha;
        }
        else if(spline_type == Vanderbilt)
        {
            // At r = Rmax, these boundary conditions force f(r) = g(r) (effective mass->infinity).
            A[(n-1)*n2 + (n-1)] += 0.5/Constant::Alpha;
            A[(n-1)*n2 + (n2-1)] += -0.5;
            A[(n2-1)*n2 + (n-1)] += 0.5;
            A[(n2-1)*n2 + (n2-1)] += -0.5*Constant::Alpha;
        }

        // Solve A*p[i] = energy*B*p[i]
        if(E.SolveMatrixEquation(A, b, eigenvalues, n2))
        {
            if(debug)
                *outstream << "kappa = " << kappa << std::endl;

            // Interpolate onto HF lattice and store
            unsigned int count = 0;
            unsigned int pqn = l + 1;
            i = n2/2;
            const Orbital* s;

            while((count < num_states_per_l[l]) && (i < n2))
            {
                s = core->GetState(OrbitalInfo(pqn, kappa));

                // If state is not in the open shell part, check whether it is in the core
                if(s != NULL && !core->IsOpenShellState(OrbitalInfo(pqn, kappa)))
                {   if(debug)
                    {   double diff = fabs((s->Energy() - eigenvalues[i])/s->Energy());
                        *outstream << "  " << s->Name() << " en: " << std::setprecision(8) << eigenvalues[i]
                                   << "  deltaE: " << diff << std::endl;
                    }
                }
                else
                {   Orbital* ds = new Orbital(pqn, kappa);
                    ds->SetEnergy(eigenvalues[i]);
                    ds->ReSize(lattice_size);

                    double* coeff = &A[i*n2];

                    for(j=0; j<n2; j++)
                    {
                        for(point = 0; point < lattice_size; point++)
                        {
                            ds->f[point] += *coeff * B_lattice[j].f[point];
                            ds->g[point] += *coeff * B_lattice[j].g[point];
                            ds->df[point] += *coeff * B_lattice[j].df[point] * dR_lattice[point];
                            ds->dg[point] += *coeff * B_lattice[j].dg[point] * dR_lattice[point];
                        }
                        coeff++;
                    }

                    if(fabs(ds->Norm(lattice) - 1.) > 1.e-2)
                    {   if(debug)
                            *outstream << "  SingleParticleWavefunction removed: energy = " << ds->Energy()
                                       << "  norm = " << ds->Norm(lattice) << std::endl;
                        pqn--;

                        delete ds;
                    }
                    else
                    {   // Check if state already exists
                        Orbital* existing = GetState(OrbitalInfo(pqn, kappa));
                        if(existing)
                        {   *existing = *ds;
                            delete ds;
                        }
                        else
                            AddState(ds);
                        count++;
                    }
                }
                pqn++;
                i++;
            }

            // Delete unwanted higher states
            StateSet::iterator it = AllStates.find(OrbitalInfo(pqn, kappa));
            while(it != AllStates.end())
            {
                it->second.DeleteState();
                AllStates.erase(it);

                pqn++;
                it = AllStates.find(OrbitalInfo(pqn, kappa));
            }

            // Repeat orthogonalisation and print
            StateIterator basis_it = GetStateIterator();
            basis_it.First();
            while(!basis_it.AtEnd())
            {
                if(basis_it.GetState()->Kappa() == kappa)
                {
                    if(orthogonalise_again)
                        Orthogonalise(basis_it.GetState());

                    if(debug)
                    {   const Orbital* ds = basis_it.GetState();
                        s = core->GetState(OrbitalInfo(ds->RequiredPQN(), kappa));
                        if(s)
                        {   double diff = fabs((s->Energy() - ds->Energy())/s->Energy());
                            *outstream << "  " << ds->Name() << " en: " << std::setprecision(8) << ds->Energy()
                                       << " norm: " << ds->Norm(lattice) - 1. << "  deltaE: " << diff << std::endl;
                        }
                        else
                            *outstream << "  " << ds->Name() << " en: " << std::setprecision(8) << ds->Energy()
                                       << " norm: " << ds->Norm(lattice) - 1. << std::endl;
                    }
                }
                basis_it.Next();
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

    if(debug)
        *outstream << "Basis Orthogonality test: " << TestOrthogonality() << std::endl;
}

void BSplineBasis::Update()
{
    ClearSigmas();
    CreateExcitedStates(NumStatesPerL);
}
