#include "Include.h"
#include "ConfigGenerator.h"

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

void ConfigGenerator::GenerateMultipleExcitations(ConfigList& configlist, unsigned int num_excitations)
{
    ConfigList::iterator it = configlist.begin();
    Parity  parity = it->GetParity();

    for(unsigned int i=0; i<num_excitations; i++)
    {
        GenerateExcitations(configlist);
        configlist.sort();
        configlist.unique();
    }

    // Remove states of wrong parity
    it = configlist.begin();

    while(it != configlist.end())
    {   
        if(it->GetParity() != parity)
        {   //std::cout << "                              " << it->Name() << std::endl;
            it = configlist.erase(it);
        }
        else
            it++;
    }
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

void ConfigGenerator::GenerateRelativisticConfigs(const ConfigList& nrlist, RelativisticConfigList& rlist)
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

void ConfigGenerator::GenerateProjections(RelativisticConfigList& rlist, int two_m)
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
