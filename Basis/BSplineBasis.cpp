#include "BSplineBasis.h"
#include "BSplineGrid.h"
#include "BSplineCore.h"
#include "Include.h"
#include "Spline.h"
#include "Universal/CoupledFunction.h"
#include "Universal/Constant.h"
#include "Universal/Eigensolver.h"

#include "HartreeFock/StateIntegrator.h"
#include "Basis/HFExcitedStates.h"

void BSplineBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    if(!num_states_per_l.size())
        return;

    NumStatesPerL = num_states_per_l;
    Clear();

    bool debug = DebugOptions.OutputHFExcited();

    double dr0;

    Eigensolver E;

    for(unsigned int l=0; l<num_states_per_l.size(); l++)
    {
        if(num_states_per_l[l])
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

        BSplineGrid grid(n, k, dr0, rmax);
        unsigned int size = grid.Size();
        const double* R = grid.R();
        const double* dR = grid.dR();

        double* B = new double[n*size];
        double* dB = new double[n*size];
        unsigned int i, j;
        for(i=0; i<n*size; i++)
            B[i] = dB[i] = 0.;

        // Get values of splines on grid
        int nderiv = 2;     // Number of derivatives of spline (including zeroth)
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

            double x = R[point];

            // 'left' is an array marker, so add 1 for fortran arrays
            int leftplus1 = (int)left + 1;
            bsplvd_(knots, &k, &x, &leftplus1, fspline_buf, &nderiv);

            // Transfer spline values from fspline_buf
            for(unsigned int s = 0; s < k; s++)
            {   B[(s + left + 1 - k)*size + point] = fspline_buf[s];
                dB[(s + left + 1 - k)*size + point] = fspline_buf[s + MaximumK];
            }

            point++;
            m++;
        }

        delete[] fspline_buf;

        // Wavefunctions are expanded in terms of B-splines as
        //     f = Sum_i (c_i. B_i)
        //     g = Sum_i (d_i. B_i)
        // To make integrals b(i,j) = <Bi|Bj> and A(i,j) = <Bi|H|Bj>

        unsigned int n2 = 2 * n;      // Size of the matrices A and b
        double* b = new double[n2*n2];
        double* A = new double[n2*n2];
        double* eigenvalues = new double[n2];
        double* bcoef = new double[n2];

        BSplineCore spline_core(&grid, core);
        std::vector<double> potential = spline_core.GetPotential();

        for(int kappa = - int(l+1); kappa <= int(l); kappa += 2*int(l) + 1)
        {
            if(kappa == 0)
                break;

            for(i=0; i<n2*n2; i++)
                b[i] = A[i] = 0.;

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
                        double BB = B[i*size + point] * B[j*size + point] * dR[point];
                        int_1 += BB;
                        int_V += BB * potential[point];
                        int_KappaOnR += BB * kappa/R[point];
                    }

                    b[i*n2 + j] = b[j*n2 + i] = int_1;
                    b[(n+i)*n2 + (n+j)] = b[(n+j)*n2 + (n+i)] = int_1; // Constant::AlphaSquared * int_1

                    A[i*n2 + j] += (-int_V);
                    A[(n+i)*n2 + (n+j)] += -2./Constant::AlphaSquared * int_1 - int_V; // (-2. * int_1 - Constant::AlphaSquared * int_V)
                    A[(n+i)*n2 + j] += -int_KappaOnR/Constant::Alpha; // int_KappaOnR
                    A[i*n2 + (n+j)] += -int_KappaOnR/Constant::Alpha; // int_KappaOnR

                    if(i != j)
                    {
                        A[j*n2 + i] += (-int_V);
                        A[(n+j)*n2 + (n+i)] += -2./Constant::AlphaSquared * int_1 - int_V; // (-2. * int_1 - Constant::AlphaSquared * int_V)
                        A[(n+j)*n2 + i] += -int_KappaOnR/Constant::Alpha; // int_KappaOnR
                        A[j*n2 + (n+i)] += -int_KappaOnR/Constant::Alpha; // int_KappaOnR
                    }
                }
            }

            // Do the same for d/dr and exchange potential
            for(i=0; i<n; i++)
            {
                CoupledFunction exf;
                CoupledFunction exg;
                spline_core.CalculateExchange(i, kappa, true, exf);
                spline_core.CalculateExchange(i, kappa, false, exg);

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
                        int_ddR += B[i*size + point] * dB[j*size + point] * dR[point];

                        int_Eff += B[j*size + point] * exf.f[point] * dR[point];
                        int_Efg += B[j*size + point] * exf.g[point] * dR[point];// * Constant::AlphaSquared;
                        int_Egf += B[j*size + point] * exg.f[point] * dR[point];
                        int_Egg += B[j*size + point] * exg.g[point] * dR[point];// * Constant::AlphaSquared;
                    }

                    A[i*n2 + (n+j)] += int_ddR/Constant::Alpha; // -int_ddR;
                    A[(n+i)*n2 + j] += -int_ddR/Constant::Alpha;

                    A[j*n2 + i] += (- int_Eff);
                    A[(n+j)*n2 + i] +=  int_Efg * Constant::Alpha; //(- Constant::AlphaSquared * int_Efg);
                    A[j*n2 + (n+i)] +=  int_Egf * Constant::Alpha;
                    A[(n+j)*n2 + (n+i)] +=  - int_Egg * Constant::AlphaSquared; //(- Constant::AlphaSquared * int_Egg);
                }
            }

            // Account for non-zero boundary conditions on splines
            if(kappa < 0)
                A[0] += 1./Constant::Alpha;
            else
                A[0] += 2./Constant::AlphaSquared;
            A[(n-1)*n2 + (n-1)] += 0.5/Constant::Alpha; // -0.5
            A[0*n2 + n] += 0.5/Constant::Alpha; // 0.5
            A[n*n2 + 0] += -0.5/Constant::Alpha; // 0.5
            A[(n2-1)*n2 + (n2-1)] += -0.5/Constant::Alpha; // 0.5

            A[(n-1)*n2 + (n2-1)] += -0.5/Constant::Alpha;
            A[(n2-1)*n2 + (n-1)] += 0.5/Constant::Alpha;

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
                i = n;
                const DiscreteState* s;

                while((count < num_states_per_l[l]) && (i < n2))
                {
                    s = core->GetState(StateInfo(pqn, kappa));

                    // If state is not in the open shell part, check whether it is in the core
                    if(s != NULL && !core->IsOpenShellState(StateInfo(pqn, kappa)))
                    {   if(debug)
                        {   double diff = fabs((s->Energy() - eigenvalues[i])/s->Energy());
                            *outstream << "  " << s->Name() << " en: " << std::setprecision(8) << eigenvalues[i]
                                       << "  deltaE: " << diff << std::endl;
                        }
                    }
                    else
                    {   DiscreteState* ds = new DiscreteState(lattice, pqn, kappa);
                        ds->SetEnergy(eigenvalues[i]);
                        ds->ReSize(HF_size);

                        for(unsigned int j=0; j<n2; j++)
                            bcoef[j] = A[i*n2 + j];

                        int jderiv;
                        double x;

                        for(point = 0; point < HF_size; point++)
                        {
                            x = HF_R[point];

                            jderiv = 0;
                            ds->f[point] = bvalue_(knots, bcoef, &n, &k, &x, &jderiv);
                            ds->g[point] = - bvalue_(knots, bcoef + n, &n, &k, &x, &jderiv)/Constant::Alpha;

                            jderiv = 1;
                            ds->df[point] = bvalue_(knots, bcoef, &n, &k, &x, &jderiv)*HF_dR[point];
                            ds->dg[point] = - bvalue_(knots, bcoef + n, &n, &k, &x, &jderiv)*HF_dR[point]/Constant::Alpha;
                        }

                        if(fabs(ds->Norm() - 1.) > 1.e-2)
                        {   if(debug)
                                *outstream << "  State removed: energy = " << ds->Energy()
                                           << "  norm = " << ds->Norm() << std::endl;
                            pqn--;

                            delete ds;
                        }
                        else
                        {   if(debug)
                            {   if(s)
                                {   double diff = fabs((s->Energy() - eigenvalues[i])/s->Energy());
                                    *outstream << "  " << ds->Name() << " en: " << std::setprecision(8) << ds->Energy()
                                               << " norm: " << ds->Norm() - 1. << "  deltaE: " << diff << std::endl;
                                }
                                else
                                    *outstream << "  " << ds->Name() << " en: " << std::setprecision(8) << ds->Energy()
                                               << " norm: " << ds->Norm() - 1. << std::endl;
                            }
                            AddState(ds);
                            count++;
                        }
                    }
                    pqn++;
                    i++;
                }
            }
        }   // kappa loop
        
        delete[] eigenvalues;
        delete[] bcoef;
        delete[] A;
        delete[] b;
        delete[] B;
        delete[] dB;
    }   }

    if(debug)
        *outstream << "Orthogonality test: " << TestOrthogonality() << std::endl;
}

void BSplineBasis::Update()
{
    Clear();

    SigmaMap::iterator sigma = SecondOrderSigma.begin();
    while(sigma != SecondOrderSigma.end())
    {   delete sigma->second;
        sigma++;
    }
    SecondOrderSigma.clear();

    CreateExcitedStates(NumStatesPerL);
}
