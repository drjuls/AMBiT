#include "Include.h"
#include "HartreeFockBasis.h"

void HartreeFockBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    if(!num_states_per_l.size())
        return;

    Clear();

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
                    if(!core->IsOpenShellState(OrbitalInfo(pqn, kappa)))
                        s = core->GetState(OrbitalInfo(pqn, kappa));

                    if(s == NULL)
                    {   // Check if state already exists
                        s = GetState(OrbitalInfo(pqn, kappa));
                        if(s == NULL)
                        {   Orbital* ds = new Orbital(kappa, 0., pqn);
                            unsigned int loop = core->CalculateExcitedState(ds);
                            if(loop)  // tells us whether ds is pre-existing OpenShellState
                                Orthogonalise(ds);

                            AddState(ds);
                            previous_state = ds;
                        }
                        else
                        {   Orbital* ds = GetState(OrbitalInfo(pqn, kappa));
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
                    Orbital* ds = GetState(OrbitalInfo(pqn, kappa));
                    if(ds == NULL)
                    {   ds = new Orbital(kappa, 0., pqn);
                        ds->SetNu(previous_state->GetNu() + 1.);
                        AddState(ds);
                    }

                    unsigned int loop = core->CalculateExcitedState(ds);
                    if(loop)  // tells us whether ds is pre-existing OpenShellState
                        Orthogonalise(ds);

                    if(DebugOptions.OutputHFExcited())
                        *outstream << "  " << ds->Name() << " en:   " << ds->GetEnergy() << "  size:  " << ds->Size() << std::endl;

                    previous_state = ds;
                    count++;
                    pqn++;
                }
            }
        }
    }

    if(DebugOptions.OutputHFExcited())
        *outstream << "Basis Orthogonality test: " << TestOrthogonality() << std::endl;
}

/** Update all of the excited states because the core has changed. */
void HartreeFockBasis::Update()
{
    ClearSigmas();

    StateIterator it = GetStateIterator();
    it.First();
    while(!it.AtEnd())
    {
        Orbital* ds = it.GetState();
        core->UpdateExcitedState(ds);
        it.Next();
    }
}
