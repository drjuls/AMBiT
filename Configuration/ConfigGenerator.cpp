#include "Include.h"
#include "ConfigGenerator.h"
#include "HartreeFock/ConfigurationParser.h"
#include <numeric>

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

pConfigList ConfigGenerator::GenerateNonRelConfigurations(pHFOperator one_body, pHartreeY two_body)
{
    user_input.set_prefix("CI");
    pConfigList biglist = ParseAndGenerateNonRelConfigurations(one_body, two_body);
    user_input.set_prefix("");

    return biglist;
}

pConfigList ConfigGenerator::GenerateNonRelConfigurations(unsigned int& small_size, pHFOperator one_body, pHartreeY two_body)
{
    user_input.set_prefix("CI");
    pConfigList biglist = ParseAndGenerateNonRelConfigurations(one_body, two_body);
    small_size = biglist->size();
    if(user_input.SectionExists("CI/SmallSide"))
    {
        user_input.set_prefix("CI/SmallSide");
        pConfigList sublist = ParseAndGenerateNonRelConfigurations(one_body, two_body);

        // Join sublist and biglist by getting set difference and appending to sublist.
        small_size = sublist->size();
        pConfigList joined = std::make_shared<ConfigList>(biglist->size() + sublist->size());
        auto itsmall = std::copy(sublist->begin(), sublist->end(), joined->begin());

        auto last = std::set_difference(biglist->begin(), biglist->end(), sublist->begin(), sublist->end(), itsmall);
        joined->resize(last - joined->begin());
        biglist = joined;
    }
    user_input.set_prefix("");

    return biglist;
}

pConfigList ConfigGenerator::ParseAndGenerateNonRelConfigurations(pHFOperator one_body, pHartreeY two_body)
{
    pConfigList nrlist = std::make_shared<ConfigList>();

    int electron_excitations = 0;
    int hole_excitations = 0;
    int num_electron_excitation_inputs = user_input.vector_variable_size("ElectronExcitations");
    int num_hole_excitation_inputs = user_input.vector_variable_size("HoleExcitations");
    
    // Default if no input detected is 2 electron excitations.
    // If input is detected it can either be a number, which would set electron_excitations and use ValenceBasis to determine states to excite to,
    // otherwise if input is of the form 'X, Y, ...' where X is the pqn and Y is the string with the basis (eg. 8spdf) then the code will
    // allow X excitations to the states in Y
    if(num_electron_excitation_inputs < 1)
    {
        num_electron_excitation_inputs = 1;
        electron_excitations = 2;
    }
    else if(num_electron_excitation_inputs == 1 && user_input("ElectronExcitations", 2) >= 0)
    {
        electron_excitations = user_input("ElectronExcitations", 2);
    }
    else if(num_electron_excitation_inputs%2 == 0)  // Input should come in pairs
    {
        electron_excitations = num_electron_excitation_inputs/2;
    }
    else
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
    else if(num_hole_excitation_inputs == 1 && user_input("HoleExcitations", 0) >= 0)
    {
        hole_excitations = user_input("HoleExcitations", 0);
    }
    else if(num_hole_excitation_inputs%2 == 0)  // Input should come in pairs
    {
        hole_excitations = num_hole_excitation_inputs/2;
    }
    else
    {
        *errstream << "USAGE: CI/HoleExcitations incorrectly specified." << std::endl;
        exit(1);
    }
    
    // Get leading configurations
    leading_configs->clear();
    int numValenceElectrons = 0;
    
    int num_configs = user_input.vector_variable_size("LeadingConfigurations");
    leading_configs->reserve(num_configs);
    for(int i = 0; i < num_configs; i++)
    {
        const std::string name = user_input("LeadingConfigurations", "", i);
        NonRelConfiguration config(name);

        // Check that the configuration gels with the number of electrons
        if(i == 0)
            numValenceElectrons = config.ElectronNumber();
        else if(config.ElectronNumber() != numValenceElectrons)
        {   *errstream << "USAGE: LeadingConfiguration " << name
                       << " does not have correct number of valence electrons." << std::endl;
            exit(1);
        }
        leading_configs->push_back(config);
    }

    if(numValenceElectrons == 0)
    {   // Add vacuum state: no holes or electrons
        leading_configs->push_back(NonRelConfiguration());
    }

    SortAndUnique(leading_configs);
    nrlist->insert(nrlist->end(), leading_configs->begin(), leading_configs->end());
    SortAndUnique(nrlist);

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
                std::string CI_basis_string = user_input("ElectronExcitations", "", 2*excitation_step + 1);
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
                std::string CI_basis_string = user_input("HoleExcitations", "", 2*excitation_step + 1);
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

    // Trim configurations outside of ConfigurationAverageEnergyRange
    if(one_body && two_body && user_input.vector_variable_size("ConfigurationAverageEnergyRange") == 2)
    {
        double lower_energy = user_input("ConfigurationAverageEnergyRange", 0.0, 0);
        double upper_energy = user_input("ConfigurationAverageEnergyRange", 0.0, 1);

        auto it = nrlist->begin();
        while(it != nrlist->end())
        {
            double energy = it->CalculateConfigurationAverageEnergy(orbitals->valence, one_body, two_body);
            if(energy < lower_energy || energy > upper_energy)
                it = nrlist->erase(it);
            else
                it++;
        }
    }

    // Add extra configurations not to be considered leading configurations.
    // These are added regardless of ConfigurationAverageEnergyRange
    int num_extra_configs = user_input.vector_variable_size("ExtraConfigurations");
    nrlist->reserve(nrlist->size() + num_extra_configs);
    for(int i = 0; i < num_extra_configs; i++)
    {
        const std::string extraname = user_input("ExtraConfigurations", "", i);
        NonRelConfiguration extraconfig(extraname);
        
        if(extraconfig.ElectronNumber() != numValenceElectrons)
        {
            *errstream << "USAGE: ExtraConfiguration " << extraname
                       << " does not have correct number of valence electrons." << std::endl;
            exit(1);
        }
        nrlist->push_back(extraconfig);
    }

    SortAndUnique(nrlist);

    // Print configuration list if requested
    if(user_input.search("CI/--print-configurations"))
    {
        *outstream << "\nConfigurations:" << std::endl;
        for(auto& config: *nrlist)
        {
            *outstream << "  " << config << std::endl;
        }
    }

    // Print configuration average energies if requested
    if(one_body && two_body && user_input.search("--print-configuration-average-energy"))
    {
        *outstream << "\nConfiguration average energies:" << std::endl;
        for(auto& config: *nrlist)
        {
            *outstream << config << ": "
                       << config.CalculateConfigurationAverageEnergy(orbitals->valence, one_body, two_body)
                       << std::endl;
        }
    }

    return nrlist;
}

pRelativisticConfigList ConfigGenerator::GenerateRelativisticConfigurations(pConfigList nrlist) const
{
    pRelativisticConfigList rlist = std::make_shared<RelativisticConfigList>();
    ConfigList::iterator it = nrlist->begin();
    while(it != nrlist->end())
    {
        rlist->append(*it->GenerateRelativisticConfigs());
        it++;
    }

    rlist->sort(FewestProjectionsFirstComparator());
    rlist->unique();
    return rlist;
}

pRelativisticConfigList ConfigGenerator::GenerateRelativisticConfigurations(pConfigList nrlist, const Symmetry& sym, pAngularDataLibrary angular_library) const
{
    unsigned int end = nrlist->size();
    return GenerateRelativisticConfigurations(nrlist, end, sym, angular_library);
}

pRelativisticConfigList ConfigGenerator::GenerateRelativisticConfigurations(pConfigList nrlist, unsigned int& small_size, const Symmetry& sym, pAngularDataLibrary angular_library) const
{
    pRelativisticConfigList prlist = std::make_shared<RelativisticConfigList>();
    RelativisticConfigList& rlist = *prlist;

    auto check_symmetry = [&sym](const NonRelConfiguration& nrconfig) {
        return (nrconfig.GetParity() == sym.GetParity() && sym.GetTwoJ() <= nrconfig.GetTwiceMaxProjection());
    };

    // [0, small_size)
    rlist = std::accumulate(nrlist->begin(), std::next(nrlist->begin(), small_size), rlist,
                            [&check_symmetry](RelativisticConfigList& list, NonRelConfiguration& item) {
                                if(check_symmetry(item))
                                    list.append(*item.GenerateRelativisticConfigs());
                                return list;
                            });

    unsigned int rel_small_size = rlist.size();

    // [small_size, end)
    rlist = std::accumulate(std::next(nrlist->begin(), small_size), nrlist->end(), rlist,
                            [&check_symmetry](RelativisticConfigList& list, NonRelConfiguration& item) {
                                if(check_symmetry(item))
                                    list.append(*item.GenerateRelativisticConfigs());
                                return list;
                            });

    prlist->SetSmallSize(rel_small_size);

    if(angular_library)
        GenerateProjections(prlist, sym.GetTwoJ(), sym.GetTwoJ(), angular_library);

    return prlist;
}

pRelativisticConfigList ConfigGenerator::GenerateRelativisticConfigurations(const Symmetry& sym, pAngularDataLibrary angular_library)
{
    unsigned int small_size;
    pConfigList nrlist = GenerateNonRelConfigurations(small_size);

    return GenerateRelativisticConfigurations(nrlist, small_size, sym, angular_library);
}

void ConfigGenerator::GenerateExcitations(pConfigList configlist, const NonRelInfoSet& electron_valence, const NonRelInfoSet& hole_valence) const
{
    ConfigList old_list(*configlist);

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
                            configlist->push_back(new_config);
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
                            configlist->push_back(new_config);
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
                        configlist->push_back(new_config);
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
                    configlist->push_back(new_config);
                }
                else if((particle_it->second > 0) && (other_particle_it->second < 0))
                {
                    NonRelConfiguration new_config(*config_it);
                    new_config.RemoveSingleParticle(particle_it->first);
                    new_config.AddSingleParticle(other_particle_it->first);
                    configlist->push_back(new_config);
                }

                other_particle_it++;
            }

            particle_it++;
        }

        config_it++;
    }

    SortAndUnique(configlist);
}

void ConfigGenerator::GenerateProjections(pRelativisticConfigList rlist, const Symmetry& sym, int two_m, pAngularDataLibrary angular_library) const
{
    if(rlist == nullptr || rlist->size() == 0)
        return;

    // Get correct library
    auto it = rlist->begin();
    while(it != rlist->end())
    {
        if(it->GetProjections(angular_library, sym, two_m))
            it++;
        else
            it = rlist->erase(it);
    }

    angular_library->GenerateCSFs();

    // Write even if there are no CSFs for a given J since this is not so obvious
    angular_library->Write();

    // Remove from list if there are no CSFs for a particular RelativisticConfiguration.
    it = rlist->begin();
    while(it != rlist->end())
    {
        if(it->NumCSFs())
            it++;
        else
            it = rlist->erase(it);
    }
}
