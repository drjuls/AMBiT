#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "ConfigGenerator.h"
#include "HartreeFock/ConfigurationParser.h"

ConfigGenerator::ConfigGenerator(pOrbitalManagerConst orbitals, MultirunOptions& userInput):
    orbitals(orbitals), user_input(userInput)
{
    NonRelSet.clear();
    leading_configs = pConfigList(new ConfigList());

    pOrbitalMapConst valence = orbitals->GetOrbitalMap(OrbitalClassification::valence);
    auto it = valence->begin();
    while(it != valence->end())
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

pConfigList ConfigGenerator::GetLeadingConfigs()
{   return leading_configs;
}

pConfigListConst ConfigGenerator::GetLeadingConfigs() const
{   return leading_configs;
}

pRelativisticConfigList ConfigGenerator::GenerateRelativisticConfigurations(const Symmetry& sym, bool generate_projections)
{
    int electron_excitations = 0;
    int hole_excitations = 0;
    int num_electron_excitation_inputs = user_input.vector_variable_size("CI/ElectronExcitations");
    int num_hole_excitation_inputs = user_input.vector_variable_size("CI/HoleExcitations");

    // Default if no input detected is 2 electron excitations.
    // If input is detected it can either be a number, which would set electron_excitations and use ValenceBasis to determine states to excite to,
    // otherwise if input is of the form 'X, Y, ...' where X is the pqn and Y is the string with the basis (eg. 8spdf) then the code will
    // allow X excitations to the states in Y
    if(num_electron_excitation_inputs < 1)
    {
        num_electron_excitation_inputs = 1;
        electron_excitations = 2;
    }
    else if(num_electron_excitation_inputs == 1 && user_input("CI/ElectronExcitations", 2) >= 0)
    {
        electron_excitations = user_input("CI/ElectronExcitations", 2);
    }
    else if(num_electron_excitation_inputs%2 != 0)  // Input should come in pairs
    {
        *errstream << "USAGE: CI/ElectronExcitations incorrectly specified." << std::endl;
        exit(1);
    }

    // Holes have default 0 excitations.
    if(num_hole_excitation_inputs < 1)
    {
        num_hole_excitation_inputs = 1;
        hole_excitations = 0;
    }
    else if(num_hole_excitation_inputs == 1 && user_input("CI/HoleExcitations", 0) >= 0)
    {
        hole_excitations = user_input("CI/HoleExcitations", 0);
    }
    else if(num_hole_excitation_inputs%2 != 0)  // Input should come in pairs
    {
        *errstream << "USAGE: CI/HoleExcitations incorrectly specified." << std::endl;
        exit(1);
    }

    // Get leading configurations
    leading_configs->clear();
    int numValenceElectrons = 0;

    int num_configs = user_input.vector_variable_size("CI/LeadingConfigurations");
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

    if(numValenceElectrons == 0)
    {   // Add vacuum state: no holes or electrons
        leading_configs->add(NonRelConfiguration());
    }

    ConfigList nrlist;
    nrlist.merge(*leading_configs);

    // Total number of excitation steps
    int total_num_excitations = mmax(electron_excitations, hole_excitations);

    for(int excitation_step = 0; excitation_step < total_num_excitations; excitation_step++)
    {
        NonRelInfoSet valence_electrons;
        NonRelInfoSet valence_holes;

        // Electron subset
        if(excitation_step < electron_excitations)
        {
            for(const auto& pair: *orbitals->particle)
                valence_electrons.insert(NonRelInfo(pair.first));

            if(num_electron_excitation_inputs > 1)
            {
                std::string CI_basis_string = user_input("CI/ElectronExcitations", "", 2*excitation_step + 1);
                std::vector<int> limits = ConfigurationParser::ParseBasisSize(CI_basis_string);

                auto valence_it = valence_electrons.begin();
                while(valence_it != valence_electrons.end())
                {
                    int L = valence_it->L();

                    if(L >= limits.size())
                        valence_it = valence_electrons.erase(valence_it);
                    else if(valence_it->PQN() > limits[L])
                        valence_it = valence_electrons.erase(valence_it);
                    else
                        valence_it++;
                }
            }
        }

        // Holes subset
        if(excitation_step < hole_excitations)
        {
            for(const auto& pair: *orbitals->hole)
                valence_holes.insert(NonRelInfo(pair.first));

            if(num_hole_excitation_inputs > 1)
            {
                std::string CI_basis_string = user_input("CI/HoleExcitations", "", 2*excitation_step + 1);
                std::vector<int> limits = ConfigurationParser::ParseBasisSize(CI_basis_string);

                auto valence_it = valence_holes.begin();
                while(valence_it != valence_holes.end())
                {
                    int L = valence_it->L();

                    if(L >= limits.size())
                        valence_it = valence_holes.erase(valence_it);
                    else if(valence_it->PQN() > limits[L])
                        valence_it = valence_holes.erase(valence_it);
                    else
                        valence_it++;
                }
            }
        }

        GenerateExcitations(nrlist, valence_electrons, valence_holes);
    }

    // Add extra configurations not to be considered leading configurations.
    int num_extra_configs = user_input.vector_variable_size("CI/ExtraConfigurations");
    for(int i = 0; i < num_extra_configs; i++)
    {
        const std::string extraname = user_input("CI/ExtraConfigurations", "", i);
        NonRelConfiguration extraconfig(extraname);

        if(extraconfig.ElectronNumber() != numValenceElectrons)
        {
            *errstream << "USAGE: ExtraConfiguration " << extraname
            << " does not have correct number of valence electrons." << std::endl;
            exit(1);
        }
        nrlist.add(extraconfig);
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

    pRelativisticConfigList rlist(GenerateRelativisticConfigs(nrlist));

    if(generate_projections)
        GenerateProjections(rlist, sym.GetTwoJ());

    return rlist;
}

void ConfigGenerator::GenerateExcitations(ConfigList& configlist, const NonRelInfoSet& electron_valence, const NonRelInfoSet& hole_valence) const
{
    ConfigList old_list(configlist);

    // Go through the set of initial configurations
    ConfigList::const_iterator config_it = old_list.begin();
    while(config_it != old_list.end())
    {
        // Move electrons and holes:
        // For each single particle state in the configuration
        NonRelConfiguration::const_iterator particle_it = config_it->begin();
        while(particle_it != config_it->end())
        {
            // Electron or hole?
            if(particle_it->second > 0)
            {
                // Get another single particle state to move to
                for(const auto& electron: electron_valence)
                {
                    if(electron != particle_it->first)
                    {
                        NonRelConfiguration new_config(*config_it);
                        new_config.RemoveSingleParticle(particle_it->first);
                        if(new_config.AddSingleParticle(electron))
                            configlist.add(new_config);
                    }
                }
            }
            else
            {   // Get another single particle state to move to
                for(const auto& hole: hole_valence)
                {
                    if(hole != particle_it->first)
                    {
                        NonRelConfiguration new_config(*config_it);
                        new_config.AddSingleParticle(particle_it->first);
                        if(new_config.RemoveSingleParticle(hole))
                            configlist.add(new_config);
                    }
                }
            }

            particle_it++;
        }

        // Pair creation
        for(const auto& electron: electron_valence)
        {
            NonRelConfiguration config_with_extra_electron(*config_it);
            if(config_with_extra_electron.AddSingleParticle(electron))
            {
                for(const auto& hole: hole_valence)
                {
                    NonRelConfiguration new_config(config_with_extra_electron);
                    if(new_config.RemoveSingleParticle(hole))
                        configlist.add(new_config);
                }
            }
        }

        // Pair annihilation
        particle_it = config_it->begin();
        while(particle_it != config_it->end())
        {
            auto other_particle_it = particle_it;
            other_particle_it++;

            while(other_particle_it != config_it->end())
            {
                if((particle_it->second < 0) && (other_particle_it->second > 0))
                {
                    NonRelConfiguration new_config(*config_it);
                    new_config.AddSingleParticle(particle_it->first);
                    new_config.RemoveSingleParticle(other_particle_it->first);
                    configlist.add(new_config);
                }
                else if((particle_it->second > 0) && (other_particle_it->second < 0))
                {
                    NonRelConfiguration new_config(*config_it);
                    new_config.RemoveSingleParticle(particle_it->first);
                    new_config.AddSingleParticle(other_particle_it->first);
                    configlist.add(new_config);
                }

                other_particle_it++;
            }

            particle_it++;
        }

        config_it++;
    }

    configlist.unique();
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
    Parity parity = rlist->begin()->GetParity();

    std::string angular_directory = string_macro(ANGULAR_DATA_DIRECTORY);
    if(user_input.search("AngularDataDirectory"))
        angular_directory = user_input("AngularDataDirectory", "");

    // Map particle number to angular data library
    std::map<int, pAngularDataLibrary> angular_libraries;

    auto it = rlist->begin();
    while(it != rlist->end())
    {
        // Get correct library
        auto lib_iterator = angular_libraries.find(it->ParticleNumber());
        pAngularDataLibrary lib;

        if(lib_iterator == angular_libraries.end())
        {   lib.reset(new AngularDataLibrary(it->ParticleNumber(), Symmetry(two_m, parity), two_m, angular_directory));
            lib->Read();
            angular_libraries[it->ParticleNumber()] = lib;
        }
        else
            lib = lib_iterator->second;

        if(it->GetProjections(lib))
            it++;
        else
            it = rlist->erase(it);
    }

    for(auto lib_iterator: angular_libraries)
    {
        lib_iterator.second->GenerateCSFs();

        // Write even if there are no CSFs for a given J since this is not so obvious
        lib_iterator.second->Write();
    }

    it = rlist->begin();
    while(it != rlist->end())
    {
        if(it->NumCSFs())
            it++;
        else
            it = rlist->erase(it);
    }
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
        int num_particles = abs(num_electrons);
        int start = mmax(0, num_particles - rinfo2.MaxNumElectrons());
        int end = mmin(num_particles, rinfo1.MaxNumElectrons());

        // Holes
        if(num_electrons < 0)
        {
            std::swap(start, end);
            start = -start;
            end = -end;
        }

        // Next orbital is the same for all loops
        current_orbital++;
        NonRelConfiguration::const_iterator next_orbital = current_orbital;

        for(int i=start; i<=end; i++)
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
