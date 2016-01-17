#include "Include.h"
#include "BasisGenerator.h"
#include "HartreeFock/HartreeFocker.h"

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
                            *logstream << "  " << s_copy->Name() << " en:   " << s_copy->Energy() << "  size:  " << s_copy->size() << std::endl;
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
                        *logstream << "  " << ds->Name() << " en:   " << ds->Energy() << "  size:  " << ds->size() << std::endl;
                }

                pqn++;
            }
        }
    }

    return excited;
}
