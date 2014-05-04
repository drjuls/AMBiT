#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "ConfigGenerator.h"
#include <sstream>

ConfigGenerator::ConfigGenerator(pOrbitalMapConst coreStates, pOrbitalMapConst excitedStates, MultirunOptions& userInput):
    core(coreStates), excited(excitedStates), user_input(userInput)
{
    NonRelSet.clear();
    leading_configs = pConfigList(new ConfigList());

    auto it = excited->begin();
    while(it != excited->end())
    {
        pOrbitalConst ds = it->second;
        if(ds != NULL)
        {
            if(ds->Kappa() < 0)
                NonRelSet.insert(NonRelInfo(ds->PQN(), ds->L()));
        }
        it++;
    }
}

ConfigGenerator::~ConfigGenerator(void)
{
    Clear();
}

void ConfigGenerator::Clear()
{
    leading_configs->clear();
}

pConfigList ConfigGenerator::GetLeadingConfigs()
{   return leading_configs;
}

pConfigListConst ConfigGenerator::GetLeadingConfigs() const
{   return leading_configs;
}

//ConfigList* ConfigGenerator::GetNonRelConfigs()
//{
//    // Check to see if we haven't restored nrlist since a read.
//    if(nrlist.empty() && !rlist->empty())
//        RestoreNonRelConfigs();
//
//    return &nrlist;
//}

pRelativisticConfigList ConfigGenerator::GenerateRelativisticConfigurations(const Symmetry& sym)
{
    bool allow_different_excitations = false;
    unsigned int electron_excitations = 0;
    int num_electron_excitation_inputs = user_input.vector_variable_size("CI/ElectronExcitations");

    // Default if no input detected is 2 electron excitations.
    // If input is detected it can either be a number, which would set electron_excitations and use ValenceBasis to determine states to excite to,
    // otherwise if input is of the form 'X, Y, ...' where X is the pqn and Y is the string with the basis (eg. 8spdf) then the code will
    // allow X excitations to the states in Y
    if(num_electron_excitation_inputs < 1)
    {
        electron_excitations = 2;
    }
    else if(num_electron_excitation_inputs == 1 && ((int) user_input("CI/ElectronExcitations", 2)) >= 0)
    {
        electron_excitations = (int) user_input("CI/ElectronExcitations", 2);
    }
    else if(num_electron_excitation_inputs%2 != 0)  // Input should come in pairs
    {
        *errstream << "USAGE: CI/ElectronExcitations incorrectly specified." << std::endl;
        exit(1);
    }
    else
    {
        allow_different_excitations = true;
    }

    // TODO: Generate non-relativistic configs from file.
    // bool GenerateFromFile = false;
    ConfigList nrlist;

    int num_configs = user_input.vector_variable_size("CI/LeadingConfigurations");
    if(num_configs < 1)
    {   *errstream << "USAGE: Need CI/LeadingConfigurations (e.g. '3d7, 4s2 3d5')" << std::endl;
        exit(1);
    }

    Clear();
    int numValenceElectrons = 0;
    for(int i = 0; i < num_configs; i++)
    {
        const std::string name = user_input("CI/LeadingConfigurations", "", i);
        NonRelConfiguration config(name);

        // Check that the configuration gels with the number of electrons
        if(i == 0)
            numValenceElectrons = config.ElectronNumber();
        else if(config.ElectronNumber() != numValenceElectrons)
        {   *errstream << "USAGE: LeadingConfiguration " << name
                       << " does not have correct number of valence electrons." << std::endl;
            exit(1);
        }
        leading_configs->add(config);
    }

    // Adds extra configurations not to be considered leading configurations.
    if(user_input.vector_variable_size("CI/ExtraConfigurations") > 0)
    {
        int num_extra_configs = user_input.vector_variable_size("CI/ExtraConfigurations");
        for(int i = 0; i < num_extra_configs; i++)
        {
            const std::string extraname = user_input("CI/ExtraConfigurations", "", i);
            NonRelConfiguration extraconfig(extraname);

            if(extraconfig.ElectronNumber() != numValenceElectrons)
            {
                *errstream << "USAGE: LeadingConfiguration " << extraname
                           << " does not have correct number of valence electrons." << std::endl;
                exit(1);
            }
            nrlist.add(extraconfig);
        }
    }

//        if(GenerateFromFile)
//        {   ConfigFileGenerator* filegenerator = dynamic_cast<ConfigFileGenerator*>(generator);
//            filegenerator->SetInputFile("PercentagesIn.txt");
//            filegenerator->ReadConfigs(0.05);
//        }
    if(allow_different_excitations)
    {
        unsigned int CI_num_excitations;
//        std::vector<unsigned int> CI_electron_excitation_states;
        std::string CI_basis_string;

        NonRelInfoSet nrset;

        for(int i = 0; i < num_electron_excitation_inputs; i += 2)
        {
            CI_num_excitations = atoi(user_input("CI/ElectronExcitations", "", i).c_str());
            CI_basis_string = user_input("CI/ElectronExcitations", "", i+1);

//            if(!ParseBasisSize(CI_basis_string.c_str(), CI_electron_excitation_states))
//            {
//                *errstream << "USAGE: CI/ElectronExcitations = " << CI_basis_string << " incorrectly specified." << std::endl;
//                exit(1);
//            }

            nrset.clear();
            nrset.AddConfigs(CI_basis_string.c_str());

            // Erase core set
            for(auto si = core->begin(); si != core->end(); si++)
            {
                nrset.erase(si->first);
            }

            nrlist.merge(GenerateMultipleExcitationsFromLeadingConfigs(CI_num_excitations, &nrset));
        }
    }
    else
    {
        nrlist.merge(GenerateMultipleExcitationsFromLeadingConfigs(electron_excitations, &NonRelSet));
    }

    // Remove states of wrong parity
    ConfigList::iterator it = nrlist.begin();

    while(it != nrlist.end())
    {
        if(it->GetParity() != sym.GetParity())
        {   it = nrlist.erase(it);
        }
        else
            it++;
    }

    nrlist.unique();

    pRelativisticConfigList rlist = GenerateRelativisticConfigs(nrlist);
    GenerateProjections(rlist, sym.GetTwoJ());

//    auto rlist_it = rlist->begin();
//    while(rlist_it != rlist->end())
//    {
//        if(rlist_it->NumCSFs() == 0)
//            rlist_it = rlist->erase(rlist_it);
//        else
//            rlist_it++;
//    }

    return rlist;
}

void ConfigGenerator::GenerateMultipleExcitations(ConfigList& configlist, unsigned int num_excitations, const NonRelInfoSet* states_to_be_excited_to) const
{
    configlist.unique();

    for(unsigned int i=0; i<num_excitations; i++)
    {
        GenerateExcitations(configlist, states_to_be_excited_to);
        configlist.unique();
    }
}

ConfigList ConfigGenerator::GenerateMultipleExcitationsFromLeadingConfigs(unsigned int num_excitations, const NonRelInfoSet* states_to_be_excited_to) const
{
    // Move leading configs to configlist
    ConfigList nrlist;
    for(auto& config: *leading_configs)
        nrlist.add(config);

    // Generate excitations
    GenerateMultipleExcitations(nrlist, num_excitations, states_to_be_excited_to);

    return nrlist;
}

void ConfigGenerator::GenerateExcitations(ConfigList& configlist, const NonRelInfoSet* AllowedStateSet) const
{
    const NonRelInfoSet* allowed_state_set = AllowedStateSet;
    if(allowed_state_set == NULL)
        allowed_state_set = &NonRelSet;

    ConfigList old_list(configlist);
    // Go through the set of initial configurations
    ConfigList::const_iterator it = old_list.begin();
    while(it != old_list.end())
    {   
        // For each single particle state in the configuration
        NonRelConfiguration start(*it);

        NonRelConfiguration::const_iterator m_it = start.begin();
        while(m_it != start.end())
        {
            NonRelConfiguration other(start);
            other.RemoveSingleParticle(m_it->first);

            // Get another single particle state to move to
            std::set<NonRelInfo>::const_iterator nrit = allowed_state_set->begin();
            while(nrit != allowed_state_set->end())
            {
                // If the new state is not the same as the old one
                if(*nrit != m_it->first)
                {
                    NonRelConfiguration new_config(other);
                    if(new_config.AddSingleParticle(*nrit))
                        configlist.add(new_config);
                }
                nrit++;
            }
            m_it++;
        }
        it++;
    }
}

pRelativisticConfigList ConfigGenerator::GenerateRelativisticConfigs(const ConfigList& nrlist) const
{
    pRelativisticConfigList rlist(new RelativisticConfigList());
    ConfigList::const_iterator it = nrlist.begin();
    while(it != nrlist.end())
    {
        NonRelConfiguration config(*it);
        RelativisticConfiguration rconfig;
        SplitNonRelInfo(config, config.begin(), rconfig, rlist);

        it++;
    }

    rlist->unique();
    return rlist;
}

void ConfigGenerator::GenerateProjections(pRelativisticConfigList rlist, int two_m) const
{
    int particle_number = rlist->begin()->ElectronNumber();
    Parity parity = rlist->begin()->GetParity();

    std::string angular_directory = string_macro(ANGULAR_DATA_DIRECTORY);
    if(user_input.search("AngularDataDirectory"))
        angular_directory = user_input("AngularDataDirectory", "");

    pAngularDataLibrary angular_library(new AngularDataLibrary(particle_number, Symmetry(two_m, parity), two_m, angular_directory));

    angular_library->Read();

    auto it = rlist->begin();
    while(it != rlist->end())
    {
        if(it->GetProjections(angular_library))
            it++;
        else
            it = rlist->erase(it);
    }

    angular_library->GenerateCSFs();

    // Write even if there are no CSFs for a given J since this is not so obvious
    angular_library->Write();
}

/*
void ConfigGenerator::GenerateProjections(int two_m)
{
    RelativisticConfigList::iterator it = rlist->begin();
    while(it != rlist->end())
    {
        if(it->GetTwiceMaxProjection() < two_m)
        {   it = rlist->erase(it);
        }
        else if(!it->GenerateProjections(two_m))
        {   it = rlist->erase(it);
        }
        else
            it++;
    }

    rlist->sort(RelConfProjectionSizeRanking());

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
    it = rlist->begin();
    while(it != rlist->end())
    {
        if(!it->GenerateJCoefficients(double(two_m)/2.))
            it = rlist->erase(it);
        else
            it++;
    }
#endif

    rlist->sort(RelConfNumJStatesRanking());
}
*/

void ConfigGenerator::SplitNonRelInfo(const NonRelConfiguration& config, NonRelConfiguration::const_iterator current_orbital, RelativisticConfiguration& relconfig, pRelativisticConfigList& rlist) const
{
    if(current_orbital == config.end())
    {   rlist->add(relconfig);
        return;
    }

    const NonRelInfo& nrinfo(current_orbital->first);

    if(nrinfo.L() == 0)
    {   relconfig.insert(std::make_pair(current_orbital->first.GetFirstRelativisticInfo(), current_orbital->second));
        current_orbital++;
        SplitNonRelInfo(config, current_orbital, relconfig, rlist);
    }
    else
    {   // rinfo1 has kappa = -(L+1). rinfo2 has kappa = L.
        OrbitalInfo rinfo1 = nrinfo.GetFirstRelativisticInfo();
        OrbitalInfo rinfo2 = nrinfo.GetSecondRelativisticInfo();

        int num_electrons = current_orbital->second;
        int start = mmin(0, abs(num_electrons-rinfo2.MaxNumElectrons()));
        int end = mmin(num_electrons, rinfo1.MaxNumElectrons());

        // Next orbital is the same for all loops
        current_orbital++;
        NonRelConfiguration::const_iterator next_orbital = current_orbital;

        for(unsigned int i=start; i<=end; i++)
        {
            RelativisticConfiguration new_rconfig(relconfig);

            if(i)
                new_rconfig.insert(std::make_pair(rinfo1, i));
            if(num_electrons - i)
                new_rconfig.insert(std::make_pair(rinfo2, num_electrons-i));

            SplitNonRelInfo(config, next_orbital, new_rconfig, rlist);
        }
    }
}

/*
void ConfigGenerator::RestoreNonRelConfigs()
{
    // Generate non-relativistic configurations
    nrlist.clear();
    RelativisticConfigList::const_iterator it = rlist->begin();
    while(it != rlist->end())
    {   nrlist.push_back(Configuration(*it));
        it++;
    }
    nrlist.sort();
    nrlist.unique();
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
        size = rlist->size();
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
*/
