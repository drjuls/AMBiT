#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "ConfigGenerator.h"
#include <sstream>

ConfigGenerator::ConfigGenerator(const ExcitedStates* manager, const std::string& atom_identifier, const Symmetry& sym):
    filename(atom_identifier), symmetry(sym)
{
    SetExcitedStates(manager);
    filename = filename + "." + symmetry.GetString() + ".configs";
}

ConfigGenerator::~ConfigGenerator(void)
{
    Clear();
}

void ConfigGenerator::SetExcitedStates(const ExcitedStates* manager)
{
    states = manager;

    NonRelSet.clear();

    ConstStateIterator it = states->GetConstStateIterator();
    while(!it.AtEnd())
    {
        const Orbital* ds = it.GetState();
        if(ds != NULL)
        {
            if(ds->Kappa() < 0)
                NonRelSet.insert(NonRelInfo(ds->RequiredPQN(), ds->L()));
        }
        it.Next();
    }
}

void ConfigGenerator::Clear()
{
    NonRelSet.clear();
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
{
    // Check to see if we haven't restored nrlist since a read.
    if(nrlist.empty() && !rlist.empty())
        RestoreNonRelConfigs();

    return &nrlist;
}

RelativisticConfigList* ConfigGenerator::GetRelConfigs()
{   return &rlist;
}

const RelativisticConfigList* ConfigGenerator::GetRelConfigs() const
{   return &rlist;
}

unsigned int ConfigGenerator::GetNumJStates() const
{
    unsigned int N = 0;

    RelativisticConfigList::const_iterator it = rlist.begin();
    while(it != rlist.end())
    {   N += it->NumJStates();
        it++;
    }
    
    return N;
}

void ConfigGenerator::AddLeadingConfiguration(const Configuration& config)
{
    leading_configs.insert(config);
}

void ConfigGenerator::AddNonRelConfiguration(const Configuration& config)
{
    nrlist.push_back(config);
    nrlist.sort();
    nrlist.unique();
}

void ConfigGenerator::AddLeadingConfigurations(const std::set<Configuration> config_set)
{
    std::set<Configuration>::const_iterator it = config_set.begin();
    while(it != config_set.end())
    {   leading_configs.insert(*it);
        it++;
    }
}

void ConfigGenerator::GenerateMultipleExcitations(ConfigList& configlist, unsigned int num_excitations, const NonRelInfoSet* states_to_be_excited_to)
{
    // Check to see if we haven't restored nrlist since a read.
    if(nrlist.empty() && !rlist.empty())
        RestoreNonRelConfigs();

    configlist.sort();
    configlist.unique();

    for(unsigned int i=0; i<num_excitations; i++)
    {
        GenerateExcitations(configlist, states_to_be_excited_to);
        configlist.sort();
        configlist.unique();
    }

    // Remove states of wrong parity
    ConfigList::iterator it = configlist.begin();

    while(it != configlist.end())
    {   
        if(it->GetParity() != symmetry.GetParity())
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

void ConfigGenerator::GenerateMultipleExcitationsFromLeadingConfigs(unsigned int num_excitations, const NonRelInfoSet* states_to_be_excited_to)
{
    // Move leading configs to configlist
    ConfigList configlist;
    
    std::set<Configuration>::const_iterator it = leading_configs.begin();
    while(it != leading_configs.end())
    {   configlist.push_back(*it);
        it++;
    }
    
    // Generate excitations
    GenerateMultipleExcitations(configlist, num_excitations, states_to_be_excited_to);
}

void ConfigGenerator::GenerateExcitations(ConfigList& configlist, const NonRelInfoSet* AllowedStateSet)
{
    const NonRelInfoSet* allowed_state_set = AllowedStateSet;
    if(allowed_state_set == NULL)
    {   if(NonRelSet.empty())
            SetExcitedStates(states);
        allowed_state_set = &NonRelSet;
    }

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
            std::set<NonRelInfo>::const_iterator nrit = allowed_state_set->begin();
            while(nrit != allowed_state_set->end())
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

void ConfigGenerator::GenerateProjections()
{
    GenerateProjections(symmetry.GetTwoJ());
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

void ConfigGenerator::Write() const
{
    if(ProcessorRank == 0)
    {
        FILE* fp = fopen(filename.c_str(), "wb");

        unsigned int size;

        // Write leading configurations
        size = leading_configs.size();
        fwrite(&size, sizeof(unsigned int), 1, fp);
        
        std::set<Configuration>::const_iterator it = leading_configs.begin();
        while(it != leading_configs.end())
        {   it->Write(fp);
            it++;
        }
        
        // Write relativistic configurations
        size = rlist.size();
        fwrite(&size, sizeof(unsigned int), 1, fp);

        RelativisticConfigList::const_iterator rit = rlist.begin();
        while(rit != rlist.end())
        {   rit->Write(fp);
            rit++;
        }
        
        fclose(fp);
    }
}

bool ConfigGenerator::Read()
{
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        return false;

    unsigned int size, i;

    // Read leading configurations
    leading_configs.clear();
    fread(&size, sizeof(unsigned int), 1, fp);

    for(i = 0; i < size; i++)
    {
        Configuration config;
        config.Read(fp);
        leading_configs.insert(config);
    }

    // Read relativistic configurations
    rlist.clear();
    fread(&size, sizeof(unsigned int), 1, fp);

    for(i = 0; i < size; i++)
    {
        RelativisticConfiguration rconfig;
        rconfig.Read(fp);
        rlist.push_back(rconfig);
    }

    rlist.sort(RelConfNumJStatesRanking());

    fclose(fp);

    return true;
}

void ConfigGenerator::RestoreNonRelConfigs()
{
    // Generate non-relativistic configurations
    nrlist.clear();
    RelativisticConfigList::const_iterator it = rlist.begin();
    while(it != rlist.end())
    {   nrlist.push_back(it->GetNonRelConfiguration());
        it++;
    }
    nrlist.sort();
    nrlist.unique();
}
