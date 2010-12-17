#include "Include.h"
#include "Core.h"
#include "Universal/Constant.h"
#include "Configuration/Configuration.h"
#include "NonRelInfo.h"

void Core::BuildFirstApproximation(std::string configuration)
{
    /** A very complex function. It does the following
        1. Assign electrons to states
        2. Give initial energies to the states
        3. Create a first approximate potential
        4. Iterate to better the approximations
        5. Introduce the nuclear potential
     */

    if(!configuration.size())
    {
        configuration = GetConfigData();
        if(!configuration.size())
        {   *errstream << "Core::BuildFirstApproximation: Failed to get intitial states from config.txt" << std::endl;
            *outstream << "Core::BuildFirstApproximation: Failed to get intitial states from config.txt" << std::endl;
            exit(1);
        }
    }

    std::string closed_shell_string;
    std::string open_shell_string;

    size_t colon_pos = configuration.find(':');
    if(colon_pos == std::string::npos)
        closed_shell_string = configuration;
    else
    {   closed_shell_string = configuration.substr(0, colon_pos);
        open_shell_string = configuration.substr(colon_pos+1, configuration.size()-colon_pos-1);
    }

    Configuration closed_shell_config(closed_shell_string);
    Configuration open_shell_config(open_shell_string);

    if(closed_shell_config.NumParticles() + open_shell_config.NumParticles() != Z - Charge)
    {   *errstream << "Core::BuildFirstApproximation: Incorrect electron count in configuration." << std::endl;
        exit(1);
    }

    // Split non-relativistic configuration into relativistic states
    closed_shell_config.First();
    while(!closed_shell_config.AtEnd())
    {
        DiscreteState* s1 = NULL;
        DiscreteState* s2 = NULL;

        NonRelInfo info(closed_shell_config.GetInfo());
        double occupancy = closed_shell_config.GetOccupancy();
        int L = info.L();

        if(L == 0)
        {   s2 = new DiscreteState(info.PQN(), -1);
            s2->SetOccupancy(double(closed_shell_config.GetOccupancy()));
        }
        else
        {   // Split electrons between subshells
            s1 = new DiscreteState(info.PQN(), L);
            s2 = new DiscreteState(info.PQN(), -(L+1));

            double occ1 = double(2*L)/double(4*L+2) * occupancy;
            double occ2 = double(2*L+2)/double(4*L+2) * occupancy;
            s1->SetOccupancy(occ1);
            s2->SetOccupancy(occ2);
        }

        if(s1)
            AddState(s1);
        if(s2)
            AddState(s2);
        closed_shell_config.Next();
    }
    
    open_shell_config.First();
    while(!open_shell_config.AtEnd())
    {
        DiscreteState* s1 = NULL;
        DiscreteState* s2 = NULL;

        NonRelInfo info(open_shell_config.GetInfo());
        double occupancy = open_shell_config.GetOccupancy();
        int L = info.L();

        if(L == 0)
        {   s2 = new DiscreteState(info.PQN(), -1);
            s2->SetOccupancy(double(open_shell_config.GetOccupancy()));
        }
        else
        {   // Split electrons between subshells
            s1 = new DiscreteState(info.PQN(), L);
            s2 = new DiscreteState(info.PQN(), -(L+1));

            double occ1 = double(2*L)/double(4*L+2) * occupancy;
            double occ2 = double(2*L+2)/double(4*L+2) * occupancy;
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
        {   AddState(s1);
            SetOpenShellState(s1, s1->Occupancy());
        }
        if(s2)
        {   AddState(s2);
            SetOpenShellState(s2, s2->Occupancy());
        }
        open_shell_config.Next();
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
            int zero_difference = 0;        // Difference between required and actual number of nodes of wavefunction

            double trial_nu = s->Nu();

            do
            {   s->Clear();

                // Attempt the local exchange approximation method.
                // This is encapsulated in the UpdateHFPotential method
                // when first_build = true.
                iterations = ConvergeStateApproximation(s, false);
                if(iterations >= StateParameters::MaxHFIterations)
                {   // Evil - this should never happen given our simple model.
                    *errstream << "    BuildFirstApproximation (" << m << ") " << s->Name()
                               << ": iterations = " << iterations << std::endl;
                    if(m == 9)
                    {   *outstream << "    BuildFirstApproximation (" << m << ") " << s->Name()
                                   << ": iterations = " << iterations << std::endl;
                        *outstream << "       ...try modifying energy/wavefunction tolerances." << std::endl;
                        exit(1);
                    }
                }

                zero_difference = s->NumNodes() + s->L() + 1 - s->RequiredPQN();

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
}

std::string Core::GetConfigData()
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

        fgets(buffer, 200, fp);
    }

    std::string ret;
    if(found)
        ret = buffer;

    fclose(fp);
    return ret;
}
