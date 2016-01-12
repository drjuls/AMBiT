#include "RPASolver.h"
#include "Include.h"
#include "Basis/BSplineBasis.h"
#include <Eigen/Eigen>

void RPASolver::SolveRPACore(pCore core, pHFOperator hf, pHyperfineRPAOperator rpa, bool include_negative_basis)
{
    // Get basis for DeltaOrbitals: map from kappa to complete basis
    basis.clear();
    pLattice lattice = hf->GetLattice();
    pOPIntegrator integrator = hf->GetOPIntegrator();

    double Rmax = lattice->R(core->LargestOrbitalSize());
    BSplineBasis bspline_maker(lattice, 40, 7, Rmax);

    // Go through all deltaOrbitals and make sure there is a basis set for that kappa
    for(auto& pair: *core)
    {
        pRPAOrbital orbital = std::dynamic_pointer_cast<RPAOrbital>(pair.second);
        if(orbital)
        {
            for(auto& deltapsi: orbital->deltapsi)
            {
                int kappa = deltapsi.first;
                if(basis.count(kappa) == 0)
                {
                    pOrbitalMap spline_basis;
                    if(include_negative_basis)
                        spline_basis = bspline_maker.GenerateCompleteBasis(hf, kappa);
                    else
                        spline_basis = bspline_maker.GeneratePositiveBasis(hf, kappa);

                    basis[kappa] = spline_basis;
                }
            }
        }
    }

    // Iterate RPA orbitals. At each step:
    // 1. Get deltaOrbitals for each wavefunction.
    // 2. Update potentials.

    bool debug = DebugOptions.LogHFIterations();

    pCore next_states(core->Clone());

    double deltaE, max_deltaE;
    double max_norm;
    unsigned int loop = 0;

    rpa->SetCore(core);

    do
    {   loop++;
        max_deltaE = 0.;
        max_norm = 0.;

        if(debug)
            *logstream << "RPA Iteration: " << loop << std::endl;

        // Calculate new states
        for(auto pair: *next_states)
        {
            pRPAOrbital rpa_orbital = std::dynamic_pointer_cast<RPAOrbital>(pair.second);
            if(rpa_orbital)
            {
                if(debug)
                    *logstream << "  RPA orbital " << std::setw(4) << rpa_orbital->Name() << std::endl;

                for(auto deltapsi: rpa_orbital->deltapsi)
                {
                    pDeltaOrbital orbital = deltapsi.second;

                    double old_energy = orbital->DeltaEnergy();
                    deltaE = IterateDeltaOrbital(orbital, rpa);

                    double norm = orbital->Norm(integrator);
                    max_norm = mmax(norm, max_norm);

                    if(orbital->Kappa() == rpa_orbital->Kappa())
                    {
                        max_deltaE = mmax(fabs(deltaE), max_deltaE);
                        if(debug)
                            *logstream << "    kappa = " << std::setw(3) << orbital->Kappa()
                                       << "  DE = " << std::setprecision(12) << old_energy
                                       << "  deltaDE = " << std::setprecision(4) << deltaE
                                       << "  size: (" << orbital->size()
                                       << ") " << lattice->R(orbital->size())
                                       << "  norm = " << orbital->Norm(integrator) << std::endl;
                    }
                    else if(debug)
                    {   *logstream << "    kappa = " << std::setw(3) << orbital->Kappa()
                                   << "  size: (" << orbital->size()
                                   << ") " << lattice->R(orbital->size())
                                   << "  norm = " << orbital->Norm(integrator) << std::endl;
                    }
                }
            }
        }

        // Try to remove spin-orbit partner from basis if DeltaOrbitals are growing out of control
        if(!remove_spin_orbit_partner && (max_norm > 0.1))
        {
            remove_spin_orbit_partner = true;

            // Start again
            if(debug)
                *logstream << "\nRestarting RPA, removing spin-orbit partners." << std::endl;

            loop = 0;
            for(auto pair: *next_states)
            {
                pRPAOrbital rpa_orbital = std::dynamic_pointer_cast<RPAOrbital>(pair.second);
                if(rpa_orbital)
                {
                    for(auto deltapsi: rpa_orbital->deltapsi)
                    {
                        pDeltaOrbital orbital = deltapsi.second;
                        orbital->Clear();
                        orbital->SetDeltaEnergy(0.0);
                    }
                }
            }
        }

        // Copy new states
        core.reset(next_states->Clone());

        // Update potential
        rpa->SetCore(core);

    }while((max_deltaE > EnergyTolerance) && (loop < MaxRPAIterations));
}

double RPASolver::CalculateRPAExcited(pRPAOrbital orbital, pHyperfineRPAOperator rpa)
{
    MathConstant* math = MathConstant::Instance();
    pOPIntegrator integrator = rpa->GetOPIntegrator();

    *logstream << "RPA orbital " << std::setw(4) << orbital->Name() << std::endl;

    for(int twoj = mmax(1, orbital->TwoJ()-2); twoj <= orbital->TwoJ()+2; twoj+=2)
    {
        int kappa = math->convert_to_kappa(twoj, orbital->GetParity());

        pDeltaOrbital delta = std::make_shared<DeltaOrbital>(kappa, orbital);
        orbital->deltapsi.insert(std::make_pair(kappa, delta));

        double deltaE = IterateDeltaOrbital(delta, rpa);

        if(kappa == orbital->Kappa())
        {
            *logstream << "    kappa = " << std::setw(3) << kappa
                << "  DE = " << std::setprecision(12) << deltaE
                << "  norm = " << delta->Norm(integrator) << std::endl;
        }
        else
        {   *logstream << "    kappa = " << std::setw(3) << kappa
                << "  norm = " << delta->Norm(integrator) << std::endl;
        }
    }

    return orbital->deltapsi[orbital->Kappa()]->DeltaEnergy();
}

double RPASolver::IterateDeltaOrbital(pDeltaOrbital orbital, pHyperfineRPAOperator rpa)
{
    int kappa = orbital->Kappa();
    double start_DE = orbital->DeltaEnergy();
    double delta_DE = 0.;

    pRPAOrbital parent = orbital->GetParent();
    double parent_energy = parent->Energy();
    OrbitalInfo parent_info = OrbitalInfo(parent.get());

    pOPIntegrator integrator = rpa->GetOPIntegrator();

    // Apply (f + deltaV)|a>
    SpinorFunction X_a = rpa->ApplyTo(*parent, orbital->Kappa());

    // Apply deltaEnergy term if kappa == parent->Kappa()
    if(kappa == parent->Kappa())
    {
        // Get new deltaEnergy
        double new_DE = integrator->GetInnerProduct(*parent, X_a);
        orbital->SetDeltaEnergy(new_DE);
        delta_DE = new_DE - start_DE;
    }

    orbital->Clear();

    for(auto& basis_pair: *basis[kappa])
    {
        // Remove parent from basis if necessary since deltaOrbital must be orthogonal to parent
        if((parent_info.PQN() != basis_pair.first.PQN()) ||
           (remove_spin_orbit_partner && (parent_info.L() != basis_pair.first.L())) ||
           (!remove_spin_orbit_partner && (parent_info.Kappa() != basis_pair.first.Kappa())))
        {
            pOrbitalConst beta = basis_pair.second;

            double coeff = integrator->GetInnerProduct(*beta, X_a);
            coeff = -coeff/(beta->Energy() - parent_energy);

            (*orbital) += (*beta) * coeff;
        }
    }

    return delta_DE;
}

//double RPASolver::IterateDeltaOrbital(pDeltaOrbital orbital, pHFOperator hf)
//{
//    int kappa = orbital->Kappa();
//    double start_DE = orbital->DeltaEnergy();
//
//    pRPAOrbital parent = orbital->GetParent();
//    double parent_energy = parent->Energy();
//
//    // Remove parent from basis if necessary since deltaOrbital must be orthogonal to parent
//    pOrbitalMap trimbasis = basis[kappa];
//    if(kappa == parent->Kappa())
//    {
//        auto it = trimbasis->find(OrbitalInfo(parent.get()));
//        if(it != trimbasis->end())
//            trimbasis->erase(it);
//    }
//
//    // Get coefficients of basis (x) via equation
//    //  A * x = B
//    // where A = <beta | hf | alpha>
//    //       B = E_a
//    unsigned int n = trimbasis->size();
//    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
//    Eigen::VectorXd B = Eigen::VectorXd::Constant(n, parent_energy);
//
//    pOPIntegrator integrator = hf->GetOPIntegrator();
//
//    // Get vector of DE = < parent | hf | alpha>
//    Eigen::VectorXd DE = Eigen::VectorXd::Zero(n);
//
//    auto it_alpha = trimbasis->begin();
//    for(int a = 0; a < n; a++)
//    {
//        pDeltaOrbital alpha = std::static_pointer_cast<DeltaOrbital>(it_alpha->second);
//        alpha->SetDeltaEnergy(start_DE);
//        alpha->SetParent(parent);
//        SpinorFunction H_alpha = hf->ApplyTo(*alpha);
//
//        auto it_beta = trimbasis->begin();
//        for(int b = 0; b < n; b++)
//        {
//            A(b, a) = integrator->GetInnerProduct(*it_beta->second, H_alpha);
//            ++it_beta;
//        }
//
//        // Get change in deltaEnergy
//        if(kappa == parent->Kappa())
//        {
//            DE(a) = integrator->GetInnerProduct(*parent, H_alpha);
//        }
//
//        ++it_alpha;
//    }
//
//    // Solve A x = B
//    Eigen::VectorXd x = A.colPivHouseholderQr().solve(B);
//
//    *errstream << "\nA: \n";
//    *errstream << A << std::endl;
//    *errstream << x << std::endl;
//    *errstream << "Sum of x = " << x.sum() << std::endl;
//
//    // Regenerate orbital
//    orbital->Clear();
//    it_alpha = trimbasis->begin();
//    for(int a = 0; a < n; a++)
//    {
//        (*orbital) += (*it_alpha->second) * x(a);
//        ++it_alpha;
//    }
//
//    // Adjust deltaEnergy
//    double delta_DE = 0.0;
//    if(kappa == parent->Kappa())
//    {
//        delta_DE = DE.dot(x);
//        orbital->SetDeltaEnergy(start_DE + delta_DE);
//    }
//
//    return delta_DE;
//}
