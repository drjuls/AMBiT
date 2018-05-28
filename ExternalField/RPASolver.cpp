#include "RPASolver.h"
#include "RPAOperator.h"
#include "Include.h"
#include "Basis/BSplineBasis.h"

namespace Ambit
{
void RPASolver::SolveRPACore(pHFOperatorConst hf, pRPAOperator rpa)
{
    hf0 = hf;

    // Get basis for DeltaOrbitals: map from kappa to complete basis
    basis.clear();
    pLattice lattice = hf->GetLattice();
    pIntegrator integrator = hf->GetIntegrator();

    pCore rpa_core = rpa->GetRPACore();

    // Go through all deltaOrbitals and make sure there is a basis set for that kappa
    for(auto& pair: *rpa_core)
    {
        pRPAOrbital orbital = std::dynamic_pointer_cast<RPAOrbital>(pair.second);
        if(orbital)
        {
            for(auto& deltapsi: orbital->deltapsi)
            {
                int kappa = deltapsi.first->Kappa();
                if(basis.count(kappa) == 0)
                {
                    pOrbitalMap spline_basis;
                    if(include_dirac_sea)
                        spline_basis = basis_maker->GenerateCompleteBasis(hf, kappa);
                    else
                        spline_basis = basis_maker->GeneratePositiveBasis(hf, kappa);

                    // Remove core orbitals (these effects should cancel anyway)
                    auto it = spline_basis->begin();
                    while(it != spline_basis->end())
                    {
                        if(rpa_core->GetOccupancy(it->first))
                            it = spline_basis->erase(it);
                        else
                            ++it;
                    }

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

    bool is_static = rpa->IsStaticRPA();

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
                    pDeltaOrbital orbital = deltapsi.first;
                    double old_energy = orbital->DeltaEnergy();

                    if(is_static)
                        deltaE = IterateDeltaOrbital(orbital, rpa, TDHF_propnew);
                    else
                        deltaE = IterateDeltaOrbital(deltapsi, rpa, TDHF_propnew);

                    double norm = orbital->Norm(integrator);
                    max_norm = mmax(norm, max_norm);
                    max_deltaE = mmax(fabs(deltaE), max_deltaE);

                    if(debug)
                        *logstream << "    kappa = " << std::setw(3) << orbital->Kappa()
                                   << "  DE = " << std::setprecision(12) << old_energy
                                   << "  deltaDE = " << std::setprecision(4) << deltaE
                                   << "  size: (" << orbital->size()
                                   << ") " << lattice->R(orbital->size())
                                   << "  norm = " << orbital->Norm(integrator) << std::endl;
                }
            }
        }

        // Copy new states
        rpa_core.reset(next_states->Clone());

        // Update potential
        rpa->SetRPACore(rpa_core);

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
        int kappa = deltapsi.first->Kappa();
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
        pDeltaOrbital delta = pair.first;
        double deltaE = IterateDeltaOrbital(delta, rpa, 1.);

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

    return orbital->GetDeltaPsi(orbital->Kappa())->first->DeltaEnergy();
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
        if(orbital->GetDeltaPsi(kappa) == orbital->deltapsi.end())
        {
            if(rpa->IsStaticRPA())
                orbital->deltapsi.push_back(std::make_pair(std::make_shared<DeltaOrbital>(kappa, orbital), nullptr));
            else
                orbital->deltapsi.push_back(std::make_pair(std::make_shared<DeltaOrbital>(kappa, orbital), std::make_shared<DeltaOrbital>(kappa, orbital)));
        }
    }
}

double RPASolver::IterateDeltaOrbital(pDeltaOrbital orbital, pRPAOperatorConst rpa, double propnew) const
{
    int kappa = orbital->Kappa();
    double start_DE = orbital->DeltaEnergy();
    double delta_DE = 0.;

    pRPAOrbital parent = orbital->GetParent();
    double parent_energy = parent->Energy();
    OrbitalInfo parent_info = OrbitalInfo(parent.get());

    pIntegrator integrator = rpa->GetIntegrator();

    // Apply (f + deltaV)||a>
    SpinorFunction X_a = rpa->ReducedApplyTo(*parent, orbital->Kappa());

    double new_DE = start_DE * (1. - propnew);

    // Apply deltaEnergy term if kappa == parent->Kappa()
    if(kappa == parent->Kappa())
    {
        // Get new deltaEnergy
        new_DE += propnew * integrator->GetInnerProduct(*parent, X_a);
    }

    (*orbital) *= (1.-propnew);

    for(const auto& basis_pair: *basis.at(kappa))
    {
        // Remove parent from basis if necessary since deltaOrbital must be orthogonal to parent
        if((parent_info.PQN() != basis_pair.first.PQN()) ||
           (parent_info.Kappa() != basis_pair.first.Kappa()))
        {
            pOrbitalConst beta = basis_pair.second;

            double coeff = integrator->GetInnerProduct(*beta, X_a);
            coeff = -coeff/(beta->Energy() - parent_energy);

            (*orbital) += (*beta) * coeff * propnew;

            // Get new deltaEnergy
            if(kappa != parent->Kappa())
            {
                new_DE += propnew * beta->Energy() * coeff * coeff;
            }
        }
    }

    orbital->SetDeltaEnergy(new_DE);
    delta_DE = new_DE - start_DE;

    return delta_DE;
}

double RPASolver::IterateDeltaOrbital(std::pair<pDeltaOrbital, pDeltaOrbital>& orbitals, pRPAOperatorConst rpa, double propnew) const
{
    pDeltaOrbital alpha = orbitals.first;
    pDeltaOrbital alphaplus = orbitals.second;

    int kappa = alpha->Kappa();
    double start_DE = alpha->DeltaEnergy();
    double delta_DE = 0.;

    pRPAOrbital parent = alpha->GetParent();
    double parent_energy = parent->Energy();
    OrbitalInfo parent_info = OrbitalInfo(parent.get());

    pIntegrator integrator = rpa->GetIntegrator();

    // Apply (f + deltaV)||a>
    SpinorFunction X_a = rpa->ReducedApplyTo(*parent, kappa);
    SpinorFunction Y_a = rpa->ConjugateReducedApplyTo(*parent, kappa);

    double new_DE = start_DE * (1. - propnew);

    // Get new deltaEnergy
    if(kappa == parent->Kappa())
    {
        new_DE += propnew * integrator->GetInnerProduct(*parent, X_a);
    }

    (*alpha) *= (1.-propnew);
    (*alphaplus) *= (1.-propnew);

    double omega = rpa->GetFrequency();

    for(const auto& basis_pair: *basis.at(kappa))
    {
        // Remove parent from basis if necessary since deltaOrbital must be orthogonal to parent
        if((parent_info.PQN() != basis_pair.first.PQN()) ||
           (parent_info.Kappa() != basis_pair.first.Kappa()))
        {
            pOrbitalConst beta = basis_pair.second;

            double coeff = integrator->GetInnerProduct(*beta, X_a);
            coeff = -coeff/(beta->Energy() - parent_energy - omega);
            (*alpha) += (*beta) * coeff * propnew;

            double coeffplus = integrator->GetInnerProduct(*beta, Y_a);
            coeffplus = -coeffplus/(beta->Energy() - parent_energy + omega);
            (*alphaplus) += (*beta) * coeffplus * propnew;

            // Get new deltaEnergy
            if(kappa != parent->Kappa())
            {
                new_DE += propnew * beta->Energy() * coeff * coeff;
            }
        }
    }

    alpha->SetDeltaEnergy(new_DE);
    delta_DE = new_DE - start_DE;

    return delta_DE;
}
}
