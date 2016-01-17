#include "RPASolver.h"
#include "Include.h"
#include "Basis/BSplineBasis.h"

void RPASolver::SolveRPACore(pHFOperator hf, pRPAOperator rpa)
{
    hf0 = hf;

    // Get rpa core: replace orbitals with RPAOrbitals
    pCore rpa_core = std::make_shared<Core>(hf->GetLattice());
    for(const auto& orb: *hf->GetCore())
    {
        pRPAOrbital rpa_orb = std::make_shared<RPAOrbital>(*orb.second);
        CreateDeltaOrbitals(rpa_orb, rpa);

        rpa_core->AddState(rpa_orb);
    }
    rpa_core->SetOccupancies(hf->GetCore()->GetOccupancies());

    // Get basis for DeltaOrbitals: map from kappa to complete basis
    basis.clear();
    pLattice lattice = hf->GetLattice();
    pIntegrator integrator = hf->GetIntegrator();

    // Go through all deltaOrbitals and make sure there is a basis set for that kappa
    for(auto& pair: *rpa_core)
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
                    if(include_dirac_sea)
                        spline_basis = basis_maker->GenerateCompleteBasis(hf, kappa);
                    else
                        spline_basis = basis_maker->GeneratePositiveBasis(hf, kappa);

                    basis[kappa] = spline_basis;
                }
            }
        }
    }

    // Iterate RPA orbitals. At each step:
    // 1. Get deltaOrbitals for each wavefunction.
    // 2. Update potentials.

    bool debug = DebugOptions.LogHFIterations();

    pCore next_states(rpa_core->Clone());

    double deltaE, max_deltaE;
    double max_norm;
    unsigned int loop = 0;

    rpa->SetCore(rpa_core);

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
        rpa_core.reset(next_states->Clone());

        // Update potential
        rpa->SetCore(rpa_core);

    }while((max_deltaE > EnergyTolerance) && (loop < MaxRPAIterations));
}

double RPASolver::CalculateRPAExcited(pRPAOrbital orbital, pRPAOperatorConst rpa)
{
    pIntegrator integrator = rpa->GetIntegrator();

    bool debug = DebugOptions.OutputHFExcited();
    if(debug)
        *logstream << "RPA orbital " << std::setw(4) << orbital->Name() << std::endl;

    CreateDeltaOrbitals(orbital, rpa);

    // Go through all deltaOrbitals and make sure there is a basis set for that kappa
    for(auto& deltapsi: orbital->deltapsi)
    {
        int kappa = deltapsi.first;
        if(basis.count(kappa) == 0)
        {
            pOrbitalMap spline_basis;
            if(include_dirac_sea)
                spline_basis = basis_maker->GenerateCompleteBasis(hf0, kappa);
            else
                spline_basis = basis_maker->GeneratePositiveBasis(hf0, kappa);

            basis[kappa] = spline_basis;
        }
    }

    for(auto& pair: orbital->deltapsi)
    {
        pDeltaOrbital delta = pair.second;
        double deltaE = IterateDeltaOrbital(delta, rpa);

        if(debug)
        {
            if(delta->Kappa() == orbital->Kappa())
            {
                *logstream << "    kappa = " << std::setw(3) << delta->Kappa()
                           << "  DE = " << std::setprecision(12) << deltaE
                           << "  norm = " << delta->Norm(integrator) << std::endl;
            }
            else
            {   *logstream << "    kappa = " << std::setw(3) << delta->Kappa()
                           << "  norm = " << delta->Norm(integrator) << std::endl;
            }
        }
    }

    return orbital->deltapsi[orbital->Kappa()]->DeltaEnergy();
}

void RPASolver::CreateDeltaOrbitals(pRPAOrbital orbital, pRPAOperatorConst rpa) const
{
    MathConstant* math = MathConstant::Instance();
    int twoK = 2 * rpa->GetK();
    Parity Pdelta = orbital->GetParity() * rpa->GetParity();

    orbital->deltapsi.clear();

    for(int twoj = mmax(1, orbital->TwoJ() - twoK); twoj <= orbital->TwoJ() + twoK; twoj+=2)
    {
        int kappa = math->convert_to_kappa(twoj, Pdelta);
        if(orbital->deltapsi.find(kappa) == orbital->deltapsi.end())
            orbital->deltapsi[kappa] = std::make_shared<DeltaOrbital>(kappa, orbital);
    }
}

double RPASolver::IterateDeltaOrbital(pDeltaOrbital orbital, pRPAOperatorConst rpa) const
{
    int kappa = orbital->Kappa();
    double start_DE = orbital->DeltaEnergy();
    double delta_DE = 0.;

    pRPAOrbital parent = orbital->GetParent();
    double parent_energy = parent->Energy();
    OrbitalInfo parent_info = OrbitalInfo(parent.get());

    pIntegrator integrator = rpa->GetIntegrator();

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

    for(const auto& basis_pair: *basis.at(kappa))
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
