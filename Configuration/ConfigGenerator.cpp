#include "Include.h"
#include "ConfigGenerator.h"

#ifdef _MPI
#include <mpi.h>
#endif

void ConfigGenerator::SetExcitedStates(const ExcitedStates* manager)
{
    states = manager;

    NonRelSet.clear();
    RelativisticSet.clear();

    ConstStateIterator it = states->GetConstStateIterator();
    while(!it.AtEnd())
    {
        const DiscreteState* ds = it.GetState();
        if(ds != NULL)
        {
            RelativisticSet.insert(StateInfo(ds->RequiredPQN(), ds->Kappa()));
            if(ds->Kappa() < 0)
                NonRelSet.insert(NonRelInfo(ds->RequiredPQN(), ds->L()));
        }
        it.Next();
    }

    std::set<StateInfo>::const_iterator rinfo = RelativisticSet.begin();
    unsigned int index = 0;
    while(rinfo != RelativisticSet.end())
    {
        for(int two_m = -int(rinfo->TwoJ()); two_m <= int(rinfo->TwoJ()); two_m+=2)
        {
            ElectronSet[ElectronInfo(rinfo->PQN(), rinfo->Kappa(), two_m)] = index;
            index++;
        }

        rinfo++;
    }
}

void ConfigGenerator::ClearConfigLists()
{
    leading_configs.clear();
    nrlist.clear();
    rlist.clear();
}

std::set<Configuration>* ConfigGenerator::GetLeadingConfigs()
{   return &leading_configs;
}
const std::set<Configuration>* ConfigGenerator::GetLeadingConfigs() const
{   return &leading_configs;
}

ConfigList* ConfigGenerator::GetNonRelConfigs()
{   return &nrlist;
}
const ConfigList* ConfigGenerator::GetNonRelConfigs() const
{   return &nrlist;
}

RelativisticConfigList* ConfigGenerator::GetRelConfigs()
{   return &rlist;
}
const RelativisticConfigList* ConfigGenerator::GetRelConfigs() const
{   return &rlist;
}

void ConfigGenerator::AddLeadingConfiguration(const Configuration& config)
{
    leading_configs.insert(config);
}

void ConfigGenerator::AddLeadingConfigurations(const std::set<Configuration> config_set)
{
    std::set<Configuration>::const_iterator it = config_set.begin();
    while(it != config_set.end())
    {   leading_configs.insert(*it);
        it++;
    }
}

void ConfigGenerator::GenerateMultipleExcitations(ConfigList& configlist, unsigned int num_excitations, Parity parity)
{
    configlist.sort();
    configlist.unique();

    for(unsigned int i=0; i<num_excitations; i++)
    {
        GenerateExcitations(configlist);
        configlist.sort();
        configlist.unique();
    }

    // Remove states of wrong parity
    ConfigList::iterator it = configlist.begin();

    while(it != configlist.end())
    {   
        if(it->GetParity() != parity)
        {   //std::cout << "                              " << it->Name() << std::endl;
            it = configlist.erase(it);
        }
        else
            it++;
    }
    
    // Append new configurations to nrlist.
    nrlist.insert(nrlist.end(), configlist.begin(), configlist.end());
    nrlist.sort();
    nrlist.unique();
}

void ConfigGenerator::GenerateMultipleExcitationsFromLeadingConfigs(unsigned int num_excitations, Parity parity)
{
    // Move leading configs to configlist
    ConfigList configlist;
    
    std::set<Configuration>::const_iterator it = leading_configs.begin();
    while(it != leading_configs.end())
    {   configlist.push_back(*it);
        it++;
    }
    
    // Generate excitations
    GenerateMultipleExcitations(configlist, num_excitations, parity);
}

void ConfigGenerator::GenerateExcitations(ConfigList& configlist)
{
    ConfigList old_list(configlist);
    // Go through the set of initial configurations
    ConfigList::const_iterator it = old_list.begin();
    while(it != old_list.end())
    {   
        // For each single particle state in the configuration
        Configuration start(*it);
        start.First();
        while(!start.AtEnd())
        {
            Configuration other(start);
            other.RemoveSingleParticle(start.GetInfo());

            // Get another single particle state to move to
            std::set<NonRelInfo>::const_iterator nrit = NonRelSet.begin();
            while(nrit != NonRelSet.end())
            {
                // If the new state is not the same as the old one
                if(*nrit != start.GetInfo())
                {
                    Configuration new_config(other);
                    if(new_config.AddSingleParticle(*nrit))
                        configlist.push_back(new_config);
                }
                nrit++;
            }
            start.Next();
        }
        it++;
    }
}

void ConfigGenerator::GenerateRelativisticConfigs()
{
    ConfigList::const_iterator it = nrlist.begin();
    while(it != nrlist.end())
    {
        Configuration config(*it);
        config.First();
        SplitNonRelInfo(config, rlist);

        it++;
    }
    rlist.sort();
    rlist.unique();
}

void ConfigGenerator::GenerateProjections(int two_m)
{
    RelativisticConfigList::iterator it = rlist.begin();
    while(it != rlist.end())
    {
        if(it->GetTwiceMaxProjection() < two_m)
        {   it = rlist.erase(it);
        }
        else if(!it->GenerateProjections(two_m))
        {   it = rlist.erase(it);
        }
        else
            it++;
    }

    rlist.sort(RelConfProjectionSizeRanking());

#ifdef _MPI
    const unsigned int SHARING_SIZE_LIM = 200;
    // Processors take turn getting tasks, all the way down the list.
    it = rlist.begin();
    unsigned int count = 0;
    while(it != rlist.end() && (it->GetProjections().size() > SHARING_SIZE_LIM))
    {
        if(count == ProcessorRank)
            it->GenerateJCoefficients(double(two_m)/2.);

        it++;
        count++;
        if(count == NumProcessors)
            count = 0;
    }

    // Share information
    MPI::Intracomm& comm_world = MPI::COMM_WORLD;

    unsigned int num_jstates;
    double* j_coefficients;
    unsigned int coeff_size;

    it = rlist.begin();
    count = 0;
    while(it != rlist.end() && (it->GetProjections().size() > SHARING_SIZE_LIM))
    {
        if(count == ProcessorRank)
        {   // This processor calculated the jstates for this config
            num_jstates = it->NumJStates();
            comm_world.Bcast(&num_jstates, 1, MPI::UNSIGNED, count);
            
            if(num_jstates)
            {   j_coefficients = it->GetJCoefficients();
                coeff_size = num_jstates * it->GetProjections().size();
                comm_world.Bcast(j_coefficients, coeff_size, MPI::DOUBLE, count);
                it++;
            }
            else
                it = rlist.erase(it);
        }
        else
        {   // This processor needs to get the jstates
            comm_world.Bcast(&num_jstates, 1, MPI::UNSIGNED, count);

            if(num_jstates)
            {   coeff_size = num_jstates * it->GetProjections().size();
                j_coefficients = new double[coeff_size];
                comm_world.Bcast(j_coefficients, coeff_size, MPI::DOUBLE, count);
                it->SetJCoefficients(num_jstates, j_coefficients);
                it++;
            }
            else
                it = rlist.erase(it);
        }

        count++;
        if(count == NumProcessors)
            count = 0;
    }

    // Calculate the ones with size <= SHARING_SIZE_LIM
    while(it != rlist.end())
    {
        if(!it->GenerateJCoefficients(double(two_m)/2.))
            it = rlist.erase(it);
        else
            it++;
    }

#else
    it = rlist.begin();
    while(it != rlist.end())
    {
        if(!it->GenerateJCoefficients(double(two_m)/2.))
            it = rlist.erase(it);
        else
            it++;
    }
#endif

    rlist.sort(RelConfNumJStatesRanking());
}

void ConfigGenerator::SplitNonRelInfo(Configuration config, RelativisticConfigList& rlist)
{
    if(config.AtEnd())
    {   rlist.push_back(config);
        //std::cout << RelativisticConfiguration(config).Name() << std::endl;
        return;
    }

    NonRelInfo nrinfo(config.GetInfo());

    if(nrinfo.L() == 0)
    {   config.Next();
        SplitNonRelInfo(config, rlist);
    }
    else
    {   // rinfo1 has kappa = -(L+1). rinfo2 has kappa = L.
        StateInfo rinfo1 = nrinfo.GetFirstRelativisticInfo();
        StateInfo rinfo2 = nrinfo.GetSecondRelativisticInfo();

        unsigned int num_electrons = config.GetOccupancy();
        unsigned int start = mmin((unsigned int)(0), num_electrons-rinfo2.MaxNumElectrons());
        unsigned int end = mmin(num_electrons, rinfo1.MaxNumElectrons());
        for(unsigned int i=start; i<=end; i++)
        {
            RelativisticConfiguration rconfig(config);
            rconfig.SetOccupancy(rinfo1, i);
            rconfig.SetOccupancy(rinfo2, num_electrons - i);

            if(i)
                rconfig.SetIterator(rinfo1);
            else
                rconfig.SetIterator(rinfo2);

            rconfig.Next();

            SplitNonRelInfo(rconfig, rlist);
        }
    }
}
