#include "Include.h"
#include "BasisGenerator.h"
#include "HartreeFock/HartreeFocker.h"

namespace Ambit
{
pOrbitalMap BasisGenerator::GenerateHFExcited(const std::vector<int>& max_pqn)
{
    pOrbitalMap excited(new OrbitalMap(lattice));

    if(!max_pqn.size())
        return excited;

    bool debug = DebugOptions.OutputHFExcited();

    // Create Hartree-Fock solver; define integrators.
    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    HartreeFocker HF_Solver(ode_solver);

    for(int l = 0; l < max_pqn.size(); l++)
    {
        if(!max_pqn[l])
            continue;

        for(int kappa = - (l+1); kappa <= l; kappa += 2*l + 1)
        {
            if(kappa == 0)
                break;

            unsigned int pqn = l + 1;
            double nu = 0.;

            while(pqn <= max_pqn[l])
            {
                // Get first state by HF iteration
                pOrbitalConst s = open_core->GetState(OrbitalInfo(pqn, kappa));
                if(s)
                {   if(!closed_core->GetOccupancy(OrbitalInfo(pqn,kappa)))
                    {
                        pOrbital s_copy(new Orbital(s));
                        nu = s->Nu();
                        *s_copy = *s;
                        excited->AddState(s_copy);

                        if(debug)
                            *logstream << "  " << s_copy->Name() << " en:   " << std::setprecision(12) << s_copy->Energy() << "  size:  " << s_copy->size() << std::endl;
                    }
                }
                else
                {
                    pOrbital ds = pOrbital(new Orbital(kappa, pqn));
                    if(nu)
                        ds->SetNu(nu + 1./hf->GetCharge());
                    HF_Solver.CalculateExcitedState(ds, hf);
                    nu = ds->Nu();
                    excited->AddState(ds);

                    if(debug)
                        *logstream << "  " << ds->Name() << " en:   " << std::setprecision(12) << ds->Energy() << "  size:  " << ds->size() << std::endl;
                }

                pqn++;
            }
        }
    }

    return excited;
}

void BasisGenerator::UpdateHFOrbitals(const std::vector<int>& max_pqn, pOrbitalMap existing)
{
    // Create Hartree-Fock solver; define integrators.
    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    HartreeFocker hartree_focker(ode_solver);

    bool debug = DebugOptions.OutputHFExcited();

    for(int l = 0; l < max_pqn.size(); l++)
    {
        if(!max_pqn[l])
            continue;

        for(int kappa = - (l+1); kappa <= l; kappa += 2*l + 1)
        {
            if(kappa == 0)
                break;

            unsigned int pqn = l + 1;
            double nu = 0.;

            while(pqn <= max_pqn[l])
            {
                pOrbital s = existing->GetState(OrbitalInfo(pqn, kappa));
                if(s)
                {
                    pOrbital ds = s->Clone();
                    std::pair<bool, double> success = hartree_focker.ConvergeOrbitalAndExchange(ds, hf);

                    if(debug)
                        *logstream << "  " << ds->Name() << " E = " << std::setprecision(12) << ds->Energy()
                                   << "  dE = " << success.second << "  size = " << ds->size() << std::endl;

                    // Copy into the existing basis if successful, otherwise don't
                    if(success.first)
                        *s = *ds;
                    else
                        *errstream << "    BasisGenerator::UpdateHFOrbitals() did not update " << s->Name() << " orbital.\n";
                }

                pqn++;
            }
        }
    }
}
}
