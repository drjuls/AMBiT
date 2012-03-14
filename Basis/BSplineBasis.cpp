#include "BSplineBasis.h"
#include "BSplineGrid.h"
#include "BSplineCore.h"
#include "Include.h"
#include "Spline.h"
#include "Universal/Interpolator.h"
#include "Universal/CoupledFunction.h"
#include "Universal/Constant.h"
#include "Universal/Eigensolver.h"

#include "HartreeFock/StateIntegrator.h"

void BSplineBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    if(!num_states_per_l.size())
        return;

    NumStatesPerL = num_states_per_l;
    bool debug = DebugOptions.OutputHFExcited();

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

        BSplineGrid grid(num_splines, k, dr0, rmax);
        unsigned int size = grid.Size();
        const double* R_grid = grid.R();
        const double* dR_grid = grid.dR();

        // Get values of splines on grid
        double* B_gauss = new double[num_splines*size];
        double* dB_gauss = new double[num_splines*size];
        memset(B_gauss, 0, num_splines*size*sizeof(double));
        memset(dB_gauss, 0, num_splines*size*sizeof(double));

        int nderiv = 3;     // Number of derivatives of spline (including zeroth)
        double* fspline_buf = new double[nderiv * MaximumK];

        const double* knots = grid.GetSplineKnots();

        unsigned int point = 0;
        unsigned int m = 0;
        unsigned int left = k-1;

        while(point < size)
        {
            if(m == k)
            {   left++;
                m = 0;
            }

            double x = R_grid[point];

            // 'left' is an array marker, so add 1 for fortran arrays
            int leftplus1 = (int)left + 1;
            bsplvd_(knots, &k, &x, &leftplus1, fspline_buf, &nderiv);

            // Transfer spline values from fspline_buf
            for(unsigned int s = 0; s < k; s++)
            {   B_gauss[(s + left + 1 - k)*size + point] = fspline_buf[s];
                dB_gauss[(s + left + 1 - k)*size + point] = fspline_buf[s + MaximumK];
            }

            point++;
            m++;
        }

        // Get values of splines on HF lattice
        unsigned int lattice_size = core->GetLattice()->real_to_lattice(rmax);
        const double* R_lattice = core->GetLattice()->R();
        const double* dR_lattice = core->GetLattice()->dR();

        std::vector<BSpline> B_lattice(num_splines*2, BSpline(kappa, lattice_size));

        point = 0;//core->GetLattice()->real_to_lattice(dr0);
        left = k-1;

        std::vector<double> Potential(core->GetHFPotential());

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
            }

            point++;
        }

        delete[] fspline_buf;
        unsigned int i = 0;

        if(spline_type == Reno)
        {
            // Drop the splines at zero and Rmax to satisfy boundary conditions
            std::vector<BSpline>::iterator it;
            it = B_lattice.begin();
            B_lattice.erase(it, it + abs(kappa));
            it = B_lattice.begin() + n;
            B_lattice.erase(it, it + abs(kappa) + 1);        
            B_lattice.erase(B_lattice.end() - 1);

            // Renormalise splines
//            it = B_lattice.begin();
//            while(it != B_lattice.end())
//            {   it->ReNormalise(core->GetLattice());
//                it++;
//            }
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

        BSplineCore spline_core(&grid, core);
        std::vector<double> potential_grid = spline_core.GetPotential();

        unsigned int j;
        memset(b, 0, n2*n2*sizeof(double));
        memset(A, 0, n2*n2*sizeof(double));

//        for(i=0; i<n2; i++)
//        {
//            std::string filename = "spline_" + itoa(i) + ".txt";
//            B_lattice[i].Print(filename, core->GetLattice());
//        }

/*
        // Loop over Bsplines and add to A, b as needed
        for(i=0; i<n; i++)
        {
            for(j=i; j<i+k && j<n; j++)
            {
                // <Bi|Bj>; <Bi|V|Bj>; <Bi|Kappa/r|Bj>
                double int_1 = 0.;
                double int_V = 0.;
                double int_KappaOnR = 0.;

                for(point = 0; point < size; point++)
                {
                    double BB = B_gauss[i*size + point] * B_gauss[j*size + point] * dR_grid[point];
                    int_1 += BB;
                    int_V += BB * potential_grid[point];
                    int_KappaOnR += BB * kappa/R_grid[point];
                }

                b[i*n2 + j] = b[j*n2 + i] = int_1;
                b[(n+i)*n2 + (n+j)] = b[(n+j)*n2 + (n+i)] = Constant::AlphaSquared * int_1;
/*
                A[i*n2 + j] += (-int_V);
                A[(n+i)*n2 + (n+j)] += -2./Constant::AlphaSquared * int_1 - int_V; // (-2. * int_1 - Constant::AlphaSquared * int_V)
                A[(n+i)*n2 + j] += int_KappaOnR/Constant::Alpha; // int_KappaOnR
                A[i*n2 + (n+j)] += int_KappaOnR/Constant::Alpha; // int_KappaOnR

                if(i != j)
                {
                    A[j*n2 + i] += (-int_V);
                    A[(n+j)*n2 + (n+i)] += -2./Constant::AlphaSquared * int_1 - int_V; // (-2. * int_1 - Constant::AlphaSquared * int_V)
                    A[(n+j)*n2 + i] += int_KappaOnR/Constant::Alpha; // int_KappaOnR
                    A[j*n2 + (n+i)] += int_KappaOnR/Constant::Alpha; // int_KappaOnR
                }
 
            }
        }

/*
        // Do the same for d/dr and exchange potential
        for(i=0; i<n; i++)
        {
            CoupledFunction exf_lattice, exf;
            CoupledFunction exg_lattice, exg;
//                spline_core.CalculateExchange(i, kappa, true, exf);
//                spline_core.CalculateExchange(i, kappa, false, exg);

            core->CalculateExchange(B_lattice[i], kappa, exf_lattice);
            core->CalculateExchange(B_lattice[i+n], kappa, exg_lattice);

            // Interpolate exf_lattice onto spline grid
            exf.Clear();
            exf.ReSize(grid.Size());
            exg.Clear();
            exg.ReSize(grid.Size());
            Interpolator spline_interp(core->GetLattice());
            
            int order = 6;
            double deriv;
            for(unsigned int k=0; k<grid.Size(); k++)
            {   spline_interp.Interpolate(exf_lattice.f, R_grid[k], exf.f[k], deriv, order);
                spline_interp.Interpolate(exf_lattice.g, R_grid[k], exf.g[k], deriv, order);
                spline_interp.Interpolate(exg_lattice.f, R_grid[k], exg.f[k], deriv, order);
                spline_interp.Interpolate(exg_lattice.g, R_grid[k], exg.g[k], deriv, order);
            }

            for(j=0; j<n; j++)
            {
                // <Bi|d/dr|Bj>
                double int_ddR = 0.;
                double int_Eff = 0.;
                double int_Efg = 0.;
                double int_Egf = 0.;
                double int_Egg = 0.;
                
                for(point = 0; point < size; point++)
                {
                    int_ddR += B_gauss[i*size + point] * dB_gauss[j*size + point] * dR_grid[point];

                    int_Eff += B_gauss[j*size + point] * exf.f[point] * dR_grid[point];
                    int_Efg += B_gauss[j*size + point] * exf.g[point] * dR_grid[point];//Constant::AlphaSquared;
                    int_Egf += B_gauss[j*size + point] * exg.f[point] * dR_grid[point] /Constant::AlphaSquared;
                    int_Egg += B_gauss[j*size + point] * exg.g[point] * dR_grid[point] /Constant::AlphaSquared;
                }

                A[i*n2 + (n+j)] += -int_ddR/Constant::Alpha; // -int_ddR;
                A[(n+i)*n2 + j] += int_ddR/Constant::Alpha;

                A[j*n2 + i] += (- int_Eff);
                A[(n+j)*n2 + i] +=  -int_Efg * Constant::Alpha; //(- Constant::AlphaSquared * int_Efg);
                A[j*n2 + (n+i)] +=  -int_Egf * Constant::Alpha;
                A[(n+j)*n2 + (n+i)] +=  - int_Egg * Constant::AlphaSquared; //(- Constant::AlphaSquared * int_Egg);
            }
        }
*/
        StateIntegrator SI(core->GetLattice());

        for(j=0; j<n2; j++)
        {
            BSpline& Bj = B_lattice[j];

            CoupledFunction exchange;
            core->CalculateExchange(Bj, exchange);
            exchange.ReSize(Bj.Size());

            for(i=j; i<n2; i++)
            {
                BSpline& Bi = B_lattice[i];

                double En = 0.;
                for(unsigned int p=0; p<mmin(Bi.Size(), Bj.Size()); p++)
                {
                    double EQ1 = (Bj.df[p] + kappa*Bj.f[p]/R_lattice[p])
                                - (2. + Constant::AlphaSquared*Potential[p])*Bj.g[p]
                                - Constant::AlphaSquared * exchange.g[p];
                    double EQ2 = (Bj.dg[p] - kappa*Bj.g[p]/R_lattice[p]) + Potential[p]*Bj.f[p] + exchange.f[p];

                    En = En - (Bi.f[p] * EQ2 - Bi.g[p] * EQ1)*dR_lattice[p];
                }

                A[i * n2 + j] = A[j * n2 + i] = En;
                b[i * n2 + j] = b[j * n2 + i] = Bi.Overlap(Bj, core->GetLattice());
            }
        }

//        FILE* fA = fopen("Amatrix.txt", "wt");
//        FILE* fb = fopen("bmatrix.txt", "wt");
//
//        for(i = 0; i < n2; i++)
//        {
//            for(j=0; j < n2-1; j++)
//            {   fprintf(fA, "%f, ", A[i*n2+j]);
//                fprintf(fb, "%f, ", b[i*n2+j]);
//            }
//            fprintf(fA, "%f\n", A[i*n2+j]);
//            fprintf(fb, "%f\n", b[i*n2+j]);
//        }
//        fclose(fA);
//        fclose(fb);

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
            const double* HF_R = lattice->R();
            const double* HF_dR = lattice->dR();
            unsigned int HF_size = lattice->real_to_lattice(rmax) - 1;

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
                    ds->ReSize(HF_size);

                    double* coeff = &A[i*n2];

                    for(j=0; j<n2; j++)
                    {
                        for(point = 0; point < HF_size; point++)
                        {
                            ds->f[point] += *coeff * B_lattice[j].f[point];
                            ds->g[point] += *coeff * B_lattice[j].g[point];
                            ds->df[point] += *coeff * B_lattice[j].df[point] * HF_dR[point];
                            ds->dg[point] += *coeff * B_lattice[j].dg[point] * HF_dR[point];
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
                    //Orthogonalise(basis_it.GetState());

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
        delete[] B_gauss;
        delete[] dB_gauss;
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
