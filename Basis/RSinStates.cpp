#include "RSinStates.h"
#include "Include.h"

void RSinStates::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
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
                const Orbital* s;
                const Orbital* previous_state = NULL;
                while(count == 0)
                {
                    s = NULL;

                    // If state is not in the open shell part, check whether it is in the core
                    if(!core->IsOpenShellState(StateInfo(pqn, kappa)))
                        s = core->GetState(StateInfo(pqn, kappa));

                    if(s == NULL)
                    {   // Check if state already exists
                        s = GetState(StateInfo(pqn, kappa));
                        if(s == NULL)
                        {   Orbital* ds = new Orbital(pqn, kappa);
                            unsigned int loop = core->CalculateExcitedState(ds);
                            if(loop)  // tells us whether ds is pre-existing OpenShellState
                                Orthogonalise(ds);

                            AddState(ds);
                            previous_state = ds;
                        }
                        else
                        {   Orbital* ds = GetState(StateInfo(pqn, kappa));
                            unsigned int loop = core->UpdateExcitedState(ds);
                            if(loop)
                                Orthogonalise(ds);
                            previous_state = ds;
                        }
                        count++;
                    }
                    pqn++;
                }

                // Get higher states by multiplication by R or Sin(kR)
                bool TimesR = true;
                while(count < num_states_per_l[k])
                {
                    Orbital* ds = GetState(StateInfo(pqn, kappa));
                    if(ds == NULL)
                    {   ds = new Orbital(pqn, kappa);
                        AddState(ds);
                    }

                    if(TimesR)
                        MultiplyByR(previous_state, ds);
                    else
                        MultiplyBySinR(previous_state, ds);

                    TimesR = !TimesR;

                    if(DebugOptions.OutputHFExcited())
                        *outstream << "  " << ds->Name() << " en:   " << ds->Energy() << "  size:  " << ds->Size() << std::endl;

                    previous_state = ds;
                    count++;
                    pqn++;
                }

                // Delete unwanted higher states
                StateSet::iterator it = AllStates.find(StateInfo(pqn, kappa));
                while(it != AllStates.end())
                {
                    it->second.DeleteState();
                    AllStates.erase(it);

                    pqn++;
                    it = AllStates.find(StateInfo(pqn, kappa));
                }
            }
        }
    }

    if(DebugOptions.OutputHFExcited())
        *outstream << "Basis Orthogonality test: " << TestOrthogonality() << std::endl;
}
