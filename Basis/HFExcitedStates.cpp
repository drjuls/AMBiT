#include "HFExcitedStates.h"
#include "Include.h"

void HFExcitedStates::AddState(State* s)
{
    StateManager::AddState(s);

    DiscreteState* ds;
    if((ds = dynamic_cast<DiscreteState*>(s)) != NULL)
    {
        if(LastPQNForKappa.find(s->Kappa()) == LastPQNForKappa.end())
            LastPQNForKappa[s->Kappa()] = ds->RequiredPQN();
        else if(LastPQNForKappa[s->Kappa()] < ds->RequiredPQN())
            LastPQNForKappa[s->Kappa()] = ds->RequiredPQN();
    }
    else
    {   ContinuumQNs.insert(s->Nu());
    }
}

void HFExcitedStates::Write(FILE* fp) const
{
    // Count discrete states
    ConstStateIterator ex = GetConstStateIterator();
    unsigned int num_discrete = 0;
    while(!ex.AtEnd())
    {
        const DiscreteState* ds = dynamic_cast<const DiscreteState*>(ex.GetState());
        if(ds != NULL)
            num_discrete++;
        ex.Next();
    }

    // Write discrete excited states
    fwrite(&num_discrete, sizeof(unsigned int), 1, fp);
    ex.First();
    while(!ex.AtEnd())
    {
        const DiscreteState* ds = dynamic_cast<const DiscreteState*>(ex.GetState());
        if(ds != NULL)
            ds->Write(fp);
        ex.Next();
    }

    // Write continuum excited states
    unsigned int num_cont = NumStates() - num_discrete;
    fwrite(&num_cont, sizeof(unsigned int), 1, fp);
    ex.First();
    while(!ex.AtEnd())
    {
        const ContinuumState* cs = dynamic_cast<const ContinuumState*>(ex.GetState());
        if(cs != NULL)
            cs->Write(fp);
        ex.Next();
    }
}

void HFExcitedStates::Read(FILE* fp)
{
    Clear();

    unsigned int num_states, i;
    // Read excited discrete states
    fread(&num_states, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_states; i++)
    {
        DiscreteState* ds = new DiscreteState(lattice);
        ds->Read(fp);
        AddState(ds);
    }

    // Read continuum states
    fread(&num_states, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_states; i++)
    {
        ContinuumState* cs = new ContinuumState(lattice);
        cs->Read(fp);
        AddState(cs);
    }
}

void HFExcitedStates::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    if(!num_states_per_l.size())
        return;

    for(unsigned int k=0; k<num_states_per_l.size(); k++)
    {
        for(int kappa = - int(k) - 1; kappa <= int(k); kappa += 2*int(k) + 1)
        {
            if(kappa == 0)
                break;

            unsigned int count = 0;
            unsigned int pqn = k + 1;

            State* s;
            while(count < num_states_per_l[k])
            {
                s = NULL;

                // If state is not in the open shell part, check whether it is in the core
                if(!core->IsOpenShellState(StateInfo(pqn, kappa)))
                    s = core->GetState(StateInfo(pqn, kappa));

                if(s == NULL)
                {
                    // Check that it doesn't already exist
                    s = GetState(StateInfo(pqn, kappa));
                    if(s == NULL)
                    {
                        DiscreteState* ds = new DiscreteState(lattice, pqn, kappa);
                        core->CalculateExcitedState(ds);
                        AddState(ds);
                        count++;
                    }
                }
                pqn++;
            }
        }
    }
}

void HFExcitedStates::CreateContinuum(double start_nu, double end_nu, unsigned int num_nu, unsigned int num_kappa)
{
    double delta_nu = (end_nu - start_nu * double(num_nu))*2./double(num_nu * (num_nu - 1));

    for(unsigned int i=1; i<=num_nu; i++)
    {
        double nu = (double)i * start_nu + double(i*(i-1)/2) * delta_nu;

        for(unsigned int j=1; j<=num_kappa; j++)
        {
            int kappa;
            if(j%2 == 0)
                kappa = int(j)/2;
            else
                kappa = - int(j+1)/2;

            ContinuumState* s = new ContinuumState(lattice, nu, kappa);
            unsigned int it = core->CalculateExcitedState(s);
            
            if(core->GetDebugOptions().DebugHFContinuum())
                std::cout << "  " << s->Name() << " loops:   " << it << std::endl;

            AddState(s);
        }
    }
}

double HFExcitedIterator::Weight()
{
    double weight = 0.;

    DiscreteState* ds = dynamic_cast<DiscreteState*>(it->second.GetState());
    if(ds != NULL)
    {   
        std::map<int, unsigned int>::const_iterator k_it = excited_manager->LastPQNForKappa.find(ds->Kappa());
        if(k_it->second == ds->RequiredPQN())
        {
            weight = 1. + 0.5 * pow(ds->Nu(), 3.) / pow(ds->Nu() + 1./excited_manager->Charge, 2.);
        }
        else 
            weight = 1.;
    }
    else
    {   // Continuum state
        double nu = it->second.GetState()->Nu();
        std::set<double>::const_iterator current = excited_manager->ContinuumQNs.find(nu);
        std::set<double>::const_iterator next = current;
        next++;

        if(current == excited_manager->ContinuumQNs.begin())
        {   if(next != excited_manager->ContinuumQNs.end())
                weight = (*next)/2.;//1./6. * (*current) + 1./3. * (*next);
            else
                weight = 1.;
        }
        else
        {   if(next == excited_manager->ContinuumQNs.end())
            {   std::set<double>::const_iterator prev(current);
                prev--;
                weight = (*current) - (*prev)/2;//5./6. * (*current) - 1./3. * (*prev);
            }
            else
            {   std::set<double>::const_iterator prev(current);
                prev--;

                weight = ((*next) - (*prev))/2.;
            }
        }
    }
    return weight;
}

double ConstHFExcitedIterator::Weight()
{
    double weight = 0.;

    const DiscreteState* ds = dynamic_cast<const DiscreteState*>(it->second.GetState());
    if(ds != NULL)
    {
        if(excited_manager->LastPQNForKappa.find(ds->Kappa())->second == ds->RequiredPQN())
        {
            weight = 1. + 0.5 * pow(ds->Nu(), 3.) / pow(ds->Nu() + 1./excited_manager->Charge, 2.);
        }
        else 
            weight = 1.;
    }
    else
    {   // Continuum state
        double nu = it->second.GetState()->Nu();
        std::set<double>::const_iterator current = excited_manager->ContinuumQNs.find(nu);
        std::set<double>::const_iterator next = current;
        next++;

        if(current == excited_manager->ContinuumQNs.begin())
        {   if(next != excited_manager->ContinuumQNs.end())
                weight = (*next)/2.;//1./6. * (*current) + 1./3. * (*next);
            else
                weight = 1.;
        }
        else
        {   if(next == excited_manager->ContinuumQNs.end())
            {   std::set<double>::const_iterator prev(current);
                prev--;
                weight = (*current) - (*prev)/2;//5./6. * (*current) - 1./3. * (*prev);
            }
            else
            {   std::set<double>::const_iterator prev(current);
                prev--;

                weight = ((*next) - (*prev))/2.;
            }
        }
    }
    return weight;
}

void HFExcitedStates::Update()
{
    StateIterator it = GetStateIterator();
    while(!it.AtEnd())
    {
        State* old_state = it.GetState();
        State* new_state;
        if(StateInfo(old_state).Discrete())
            new_state = new DiscreteState(*dynamic_cast<DiscreteState*>(old_state));
        else
            new_state = new ContinuumState(*dynamic_cast<ContinuumState*>(old_state));

        unsigned int iterations = core->UpdateExcitedState(new_state);
        it.ReplaceState(new_state);
        it.Next();
    }

    SigmaMap::iterator sigma = SecondOrderSigma.begin();
    while(sigma != SecondOrderSigma.end())
    {   delete sigma->second;
        sigma++;
    }
    SecondOrderSigma.clear();
}
