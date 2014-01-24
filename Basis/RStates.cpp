#include "Include.h"
#include "RStates.h"
#include "Universal/MathConstant.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/StateIntegrator.h"

void RStates::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    if(!num_states_per_l.size())
        return;

    NumStatesPerL = num_states_per_l;

    for(unsigned int k=0; k<num_states_per_l.size(); k++)
    {
        if(num_states_per_l[k])
        {
            for(int kappa = - int(k) - 1; kappa <= int(k); kappa += 2*int(k) + 1)
            {
                if(kappa == 0)
                    break;

                unsigned int count = 0; // Number of states calculated so far
                unsigned int pqn = k + 1;

                // Get first state by HF iteration
                pOrbitalConst s;
                pOrbitalConst previous_state;
                while(count == 0)
                {
                    s.reset();

                    // If state is not in the open shell part, check whether it is in the core
                    if(!core->IsOpenShellState(OrbitalInfo(pqn, kappa)))
                        s = core->GetState(OrbitalInfo(pqn, kappa));

                    if(s == NULL)
                    {   // Check if state already exists
                        s = GetState(OrbitalInfo(pqn, kappa));
                        if(s == NULL)
                        {   pOrbital ds = pOrbital(new Orbital(kappa, 0., pqn));
                            unsigned int loop = core->CalculateExcitedState(ds);
                            if(loop)  // tells us whether ds is pre-existing OpenShellState
                                Orthogonalise(ds);

                            AddState(ds);
                            previous_state = ds;
                        }
                        else
                        {   pOrbital ds = GetState(OrbitalInfo(pqn, kappa));
                            unsigned int loop = core->UpdateExcitedState(ds);
                            if(loop)
                                Orthogonalise(ds);
                            previous_state = ds;
                        }
                        count++;
                    }
                    pqn++;
                }

                // Get higher states by multiplication by R
                while(count < num_states_per_l[k])
                {
                    pOrbital ds = GetState(OrbitalInfo(pqn, kappa));
                    if(ds == NULL)
                    {   ds = pOrbital(new Orbital(kappa, 0., pqn));
                        AddState(ds);
                    }

                    MultiplyByR(previous_state, ds);

                    if(DebugOptions.OutputHFExcited())
                        *outstream << "  " << ds->Name() << " en:   " << ds->GetEnergy() << "  size:  " << ds->Size() << std::endl;

                    previous_state = ds;
                    count++;
                    pqn++;
                }

                // Delete unwanted higher states
                StateSet::iterator it = AllStates.find(OrbitalInfo(pqn, kappa));
                while(it != AllStates.end())
                {
                    AllStates.erase(it);

                    pqn++;
                    it = AllStates.find(OrbitalInfo(pqn, kappa));
                }
            }
        }
    }

    if(DebugOptions.OutputHFExcited())
        *outstream << "Basis Orthogonality test: " << TestOrthogonality() << std::endl;
}

void RStates::Update()
{
    ClearSigmas();
    CreateExcitedStates(NumStatesPerL);
}
