#include "Include.h"
#include "Core.h"
#include "Universal/Constant.h"

void Core::BuildFirstApproximation()
{
    /** A very complex function. It does the following
        1. Assign electrons to states
        2. Give initial energies to the states
        3. Create a first approximate potential
        4. Iterate to better the approximations
        5. Introduce the nuclear potential
     */

    if(!GetConfigData())
    {
        *errstream << "Failed to get intitial states from config.txt" << std::endl;
        *outstream << "Failed to get intitial states from config.txt" << std::endl;

        // Get needed states
        unsigned int number_electrons = (unsigned int)(Z - Charge);
        unsigned int electron_count = 0;

        // First, get special states
        electron_count = GetSpecialStates();

        unsigned int pqn = 0;
        int kappa;
        do
        {   pqn++;
            for(int k=0; k < 2*(int)pqn-1; k++)
            {   if((double)k/2. == (double)(k/2))
                    kappa = -k/2 - 1;
                else
                    kappa = (k+1)/2;
                
                // Make sure state doesn't already exist
                if(!GetState(StateInfo(pqn, kappa)))
                {
                    DiscreteState* s = new DiscreteState(pqn, kappa);

                    if(electron_count + 2*abs(kappa) <= number_electrons)
                    {   electron_count = electron_count + 2*abs(kappa);
                    }
                    else
                    {   // Open shell
                        s->SetOccupancy(number_electrons - electron_count);
                        electron_count = number_electrons;
                    }

                    AddState(s);
                }

                if(electron_count == number_electrons)
                    break;
            }
        }
        while(electron_count < number_electrons);
    }

    // Get first approximation to state energies
    if(NumStates() > 1)
    {
        double num_states = (double)NumStates();
        StateIterator it = GetStateIterator();
        double i = 1.;
        while(!it.AtEnd())
        {   
            double ratio = (i-1.) / (num_states-1.);
            double trial_nu = (num_states - i)/((double)Z * (num_states - 1.)) + ratio;
            trial_nu = trial_nu / (-4.*ratio*ratio + 4.*ratio + 1.);
            
            (it.GetState())->SetNu(trial_nu);
            
            it.Next();
            i += 1.;
        }

        // different for some states when 56 < Z < 72
        if(((56 < Z) && (Z < 72)) || (Z == 81))
        {
            // Alter 6s state
            DiscreteState* s = GetState(StateInfo(6, -1));
            if(s != NULL)
            {   s->SetNu(1.4);
            }

            // Alter 4f state
            s = GetState(StateInfo(4, 3));
            if(s != NULL)
            {   s->SetNu(1.);
            }
            // Alter 4f* state
            s = GetState(StateInfo(4, -4));
            if(s != NULL)
            {   s->SetNu(1.);
            }
        }
    }
    else if(!Empty())
    {
        State* s = GetStateIterator().GetState();
        s->SetNu(1./Z);
    }

    // Get first approximation to potential
    double AM = 3.1,
           DD = 0.4,
           POT = 2.235;
    double HH = exp(-AM/DD);

    POT = POT * pow(Z/55., 1./3.);
    double KOP = (double)Charge;
    if(KOP == 0.)
        KOP = 1.;

    HFPotential.clear();
    HFPotential.resize(lattice->Size());

    for(unsigned int i = 0; i < HFPotential.size(); i++)
    {
        double UR = lattice->R(i)/DD;
        double ZR;
        
        if(UR <= 20.)
        {   ZR = (Z - KOP)/((1.+POT*(lattice->R(i))) * (1.+POT*(lattice->R(i))) * (HH+1.));
            ZR = ZR / (HH * exp(UR) + 1.) + KOP;
        }
        else
            ZR = KOP;

        HFPotential[i] = ZR/lattice->R(i);
    }

    LocalExchangeApproximation.clear();
    LocalExchangeApproximation.resize(HFPotential.size());

    // Iterate states, adding in local exchange approximation.
    bool debug = DebugOptions.LogFirstBuild();

    for(unsigned int m = 0; m < 10; m++)
    {
        if(debug)
            *logstream << "First Build Iteration: " << m+1 << std::endl;

        StateIterator it = GetStateIterator();
        while(!it.AtEnd())
        {
            DiscreteState* s = it.GetState();
            unsigned int iterations = 0;
            double nu_change_factor = 0.5;
            int zero_difference = 0;        // Difference between required and actual number of zeroes of wavefunction

            double trial_nu = s->Nu();

            do
            {   s->Clear();

                // Attempt the local exchange approximation method.
                // This is encapsulated in the UpdateHFPotential method
                // when first_build = true.
                iterations = ConvergeStateApproximation(s, false);
                if(iterations >= StateParameters::MaxHFIterations)
                {   // Evil - this should never happen given our simple model.
                    *errstream << "    BuildFirstApproximation: iterations = " << iterations << std::endl;
                    //PAUSE
                    //exit(1);
                }

                zero_difference = s->NumZeroes() + s->L() + 1 - s->RequiredPQN();

                if(zero_difference)
                {   
                    if(debug)
                        *logstream << "    Zero difference: " << zero_difference << "  nu = " << trial_nu << std::endl;
                    s->SetNu(trial_nu - zero_difference/abs(zero_difference) * nu_change_factor * trial_nu);
                    trial_nu = s->Nu();
                    nu_change_factor = nu_change_factor * 0.75;
                }
            }
            while(zero_difference);

            if(debug)
            {   if(DebugOptions.HartreeEnergyUnits())
                    *logstream << "  " << s->Name() << "  en: " << s->Energy() << std::endl;
                else
                    *logstream << "  " << s->Name() << " nu:   " << s->Nu() << std::endl;
            }

            it.Next();
        }

        UpdateHFPotential(0.4, true);
    }

    //Update();
}

bool Core::GetConfigData()
{
    FILE* fp = fopen("config.txt", "r");
    if(fp == NULL)
    {   *errstream << "Failed to open file \"config.txt\"" << std::endl;
        getchar();
        exit(1);
    }

    char buffer[200];

    bool found = false;
    int new_Z, new_Charge;

    while(!found && fgets(buffer, 200, fp))
    {
        // Get Z, charge
        std::stringstream ss(buffer);
        ss >> new_Z;
        ss >> new_Charge;

        if((Z == new_Z) && (Charge == new_Charge))
        {
            double radius, thickness;
            ss >> radius;
            ss >> thickness;
            NuclearRadius = radius/Constant::AtomicToFermi;
            NuclearThickness = thickness/Constant::AtomicToFermi;
            found = true;
        }
        else
            fgets(buffer, 200, fp); // skip next line
    }

    if(!found)
        return false;

    int L, pqn, occupancy;
    int electron_count = 0;

    fgets(buffer, 200, fp);
    fclose(fp);

    bool in_core = true;
    char tempbuf[20];
    unsigned int i = 0;

    // Move to first configuration
    while(!isdigit(buffer[i]) && buffer[i] != '\0')
    {   if(buffer[i] == ':')
            in_core = false;
        i++;
    }

    while(buffer[i] != '\0')
    {
        unsigned int j=0;

        // Get config
        while(isdigit(buffer[i]))
        {   tempbuf[j] = buffer[i];
            i++; j++;
        }
        tempbuf[j] = 0;
        pqn = atoi(tempbuf);

        // Get L
        for(L = 0; L < 10; L++)
        {   if(Constant::SpectroscopicNotation[L] == buffer[i])
                break;
        }
        i++;

        // Get occupancy
        j = 0;
        while(isdigit(buffer[i]))
        {   tempbuf[j] = buffer[i];
            i++; j++;
        }
        tempbuf[j] = 0;
        occupancy = atoi(tempbuf);
        electron_count += occupancy;

        DiscreteState* s1 = NULL;
        DiscreteState* s2 = NULL;

        if(L == 0)
        {   s2 = new DiscreteState(pqn, -1);
            s2->SetOccupancy(double(occupancy));
        }
        else
        {   // Split electrons between subshells
            s1 = new DiscreteState(pqn, L);
            s2 = new DiscreteState(pqn, -(L+1));

            double occ1 = double(2*L)/double(4*L+2) * double(occupancy);
            double occ2 = double(2*L+2)/double(4*L+2) * double(occupancy);
            s1->SetOccupancy(occ1);
            s2->SetOccupancy(occ2);

            // Routine for splitting electrons keeping integer occupancies
            //unsigned int occ1 = occupancy/2;
            //occ1 = mmin(occ1, unsigned int(2*L));
            //if(occ1)
            //{   s1 = new DiscreteState(pqn, L);
            //    s1->SetOccupancy(occ1);
            //}

            //s2 = new DiscreteState(pqn, -(L+1));
            //s2->SetOccupancy(mmin(occupancy-occ1, unsigned int(2*L + 2)));
        }

        if(s1)
            AddState(s1);
        if(s2)
            AddState(s2);
        
        if(!in_core)
        {   if(s1)
            {   OpenShellStates[StateInfo(s1)] = s1->Occupancy();
            }
            if(s2)
            {   OpenShellStates[StateInfo(s2)] = s2->Occupancy();
            }
        }

        // End of core?
        if(buffer[i] == ':')
        {   in_core = false;
            i++;
        }

        // Next configuration term
        while(!isdigit(buffer[i]) && buffer[i] != '\0')
            i++;
    }

    if(electron_count == Z - Charge)
        return true;
    else
    {   *errstream << "Incorrect electron count in config.txt." << std::endl;
        return false;
    }
}

unsigned int Core::GetSpecialStates()
{
    int number_electrons = int(Z - Charge);
    unsigned int number_electrons_assigned = 0;

    FILE* fp = fopen("mendel.tab", "r");
    if(fp == NULL)
    {   *errstream << "Failed to open file \"mendel.tab\"" << std::endl;
        getchar();
        exit(1);
    }
    
    char buffer[100];
    int count, kappa, zeros, occupancy;
    unsigned int L, pqn;
    while(fgets(buffer, 100, fp))
    {
        std::stringstream s(buffer);
        s >> count;

        if(count == number_electrons)
        {
            s >> kappa;
            s >> zeros;
            s >> occupancy;

            if(kappa > 0)
                L = kappa;
            else
                L = -kappa-1;

            pqn = L + zeros + 1;

            // Make state
            DiscreteState* s = new DiscreteState(pqn, kappa);
            s->SetOccupancy(occupancy);
            AddState(s);

            number_electrons_assigned += occupancy;
        }
    }
    fclose(fp);

    return number_electrons_assigned;
}
