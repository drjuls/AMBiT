#include "Include.h"
#include "RStates.h"
#include "Universal/Constant.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/StateIntegrator.h"

void RStates::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    if(!num_states_per_l.size())
        return;

    NumStatesPerL = num_states_per_l;

    // Gotta get rid of all existing states
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
                State* s;
                DiscreteState* previous_state = NULL;
                while(count == 0)
                {
                    s = NULL;

                    // If state is not in the open shell part, check whether it is in the core
                    if(!core->IsOpenShellState(StateInfo(pqn, kappa)))
                    {   s = core->GetState(StateInfo(pqn, kappa));
                        previous_state = dynamic_cast<DiscreteState*>(s);
                    }

                    if(s == NULL)
                    {
                        // Check that it doesn't already exist
                        s = GetState(StateInfo(pqn, kappa));
                        if(s == NULL)
                        {
                            DiscreteState* ds = new DiscreteState(lattice, pqn, kappa);
                            if(core->GetCharge())
                            {   unsigned int it = core->CalculateExcitedState(ds);
                                if(it)  // tells us whether ds is pre-existing OpenShellState
                                    Orthogonalise(ds);
                            }
                            else if(previous_state)
                            {   MultiplyByR(previous_state, ds);
                                if(core->GetDebugOptions().DebugHFExcited())
                                    std::cout << "  " << ds->Name() << " en:   " << ds->Energy() << "  size:  " << ds->Size() << std::endl;
                            }
                            else
                            {   std::cerr << "Cannot form HF state, kappa = " << kappa << std::endl;
                                return;
                            }
                            AddState(ds);
                            count++;
                            previous_state = ds;
                        }
                    }
                    pqn++;
                }

                // Get higher states by multiplication by R
                while(count < num_states_per_l[k])
                {
                    DiscreteState* ds = new DiscreteState(lattice, pqn, kappa);

                    MultiplyByR(previous_state, ds);

                    if(core->GetDebugOptions().DebugHFExcited())
                        std::cout << "  " << ds->Name() << " en:   " << ds->Energy() << "  size:  " << ds->Size() << std::endl;

                    AddState(ds);
                    previous_state = ds;
                    count++;
                    pqn++;
                }
            }
        }
    }
}

void RStates::Update()
{
    //ConstDiscreteStateIterator it = GetConstDiscreteStateIterator();
    //it.First();
    //unsigned int k;
    //unsigned int max_k = 0;
    //while(!it.AtEnd())
    //{   int kappa = it.GetState()->Kappa();
    //    if(kappa > 0)
    //        k = 2 * kappa - 1;
    //    else
    //        k = -2 * (kappa + 1);
    //    max_k = mmax(k, max_k);
    //    it.Next();
    //}

    //unsigned int* num_states_per_kappa = new unsigned int[max_k + 1];
    //for(k=0; k <= max_k; k++)
    //    num_states_per_kappa[k] = 0;

    //it.First();
    //while(!it.AtEnd())
    //{   int kappa = it.GetState()->Kappa();
    //    if(kappa > 0)
    //        k = 2 * kappa - 1;
    //    else
    //        k = -2 * (kappa + 1);
    //    
    //    num_states_per_kappa[k]++;
    //    it.Next();
    //}

    Clear();

    SigmaMap::iterator sigma = SecondOrderSigma.begin();
    while(sigma != SecondOrderSigma.end())
    {   delete sigma->second;
        sigma++;
    }
    SecondOrderSigma.clear();

    CreateExcitedStates(NumStatesPerL);
}
