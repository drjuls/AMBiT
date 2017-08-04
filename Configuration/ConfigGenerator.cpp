#include "Include.h"
#include "ConfigGenerator.h"
#include "HartreeFock/ConfigurationParser.h"
#include <numeric>

ConfigGenerator::ConfigGenerator(pOrbitalManagerConst orbitals, MultirunOptions& userInput):
    orbitals(orbitals), user_input(userInput)
{
    leading_configs = std::make_shared<ConfigList>();
}

pConfigList ConfigGenerator::GetLeadingConfigs()
{   return leading_configs;
}

pConfigListConst ConfigGenerator::GetLeadingConfigs() const
{   return leading_configs;
}

pRelativisticConfigList ConfigGenerator::GenerateConfigurations(pHFOperator one_body, pHartreeY two_body)
{
    // Do sublist first since leading_configs should come from CI, rather than CI/SmallSide
    pRelativisticConfigList sublist = nullptr;
    if(user_input.SectionExists("CI/SmallSide"))
    {
        user_input.set_prefix("CI/SmallSide");
        sublist = ParseAndGenerateConfigurations(one_body, two_body);
    }

    user_input.set_prefix("CI");
    pRelativisticConfigList biglist = ParseAndGenerateConfigurations(one_body, two_body);

    // Join sublist and biglist by getting set difference and appending to sublist.
    if(sublist)
    {
        pRelativisticConfigList setdiff = std::make_shared<RelativisticConfigList>(biglist->size());
        auto last = std::set_difference(biglist->begin(), biglist->end(), sublist->begin(), sublist->end(), setdiff->begin());
        setdiff->erase(last, setdiff->end());

        biglist = sublist;
        biglist->append(*setdiff);
    }

    user_input.set_prefix("");

    return biglist;
}

pRelativisticConfigList ConfigGenerator::ParseAndGenerateConfigurations(pHFOperator one_body, pHartreeY two_body)
{
    pRelativisticConfigList rlist = std::make_shared<RelativisticConfigList>();

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
    leading_configs->first.clear();
    int numValenceElectrons = 0;
    
    int num_configs = user_input.vector_variable_size("LeadingConfigurations");
    leading_configs->first.reserve(num_configs);
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
        leading_configs->first.push_back(config);
    }

    if(numValenceElectrons == 0)
    {   // Add vacuum state: no holes or electrons
        leading_configs->first.push_back(NonRelConfiguration());
    }

    leading_configs->second = leading_configs->first.size();
    SortAndUnique(leading_configs);
    rlist->append(*GenerateRelativisticConfigurations(leading_configs));

    // Total number of excitation steps
    int total_num_excitations = mmax(electron_excitations, hole_excitations);

    for(int excitation_step = 0; excitation_step < total_num_excitations; excitation_step++)
    {
        OrbitalMap valence_electrons(*orbitals->particle);
        OrbitalMap valence_holes(*orbitals->hole);

        // Electron subset
        if(excitation_step < electron_excitations)
        {
            if(num_electron_excitation_inputs > 1)
            {
                std::string CI_basis_string = user_input("ElectronExcitations", "", 2*excitation_step + 1);
                std::vector<int> limits = ConfigurationParser::ParseBasisSize(CI_basis_string);

                auto valence_it = valence_electrons.begin();
                while(valence_it != valence_electrons.end())
                {
                    int L = valence_it->first.L();
                    
                    if(L >= limits.size())
                        valence_it = valence_electrons.erase(valence_it);
                    else if(valence_it->first.PQN() > limits[L])
                        valence_it = valence_electrons.erase(valence_it);
                    else
                        valence_it++;
                }
            }
        }

        // Holes subset
        if(excitation_step < hole_excitations)
        {
            if(num_hole_excitation_inputs > 1)
            {
                std::string CI_basis_string = user_input("HoleExcitations", "", 2*excitation_step + 1);
                std::vector<int> limits = ConfigurationParser::ParseBasisSize(CI_basis_string);
                
                auto valence_it = valence_holes.begin();
                while(valence_it != valence_holes.end())
                {
                    int L = valence_it->first.L();
                    
                    if(L >= limits.size())
                        valence_it = valence_holes.erase(valence_it);
                    else if(valence_it->first.PQN() > limits[L])
                        valence_it = valence_holes.erase(valence_it);
                    else
                        valence_it++;
                }
            }
        }

        GenerateExcitations(rlist, valence_electrons, valence_holes);
    }

    // Trim configurations outside of ConfigurationAverageEnergyRange
    if(one_body && two_body && user_input.vector_variable_size("ConfigurationAverageEnergyRange") == 2)
    {
        double lower_energy = user_input("ConfigurationAverageEnergyRange", 0.0, 0);
        double upper_energy = user_input("ConfigurationAverageEnergyRange", 0.0, 1);

        auto it = rlist->begin();
        while(it != rlist->end())
        {
            double energy = it->CalculateConfigurationAverageEnergy(orbitals->valence, one_body, two_body);
            if(energy < lower_energy || energy > upper_energy)
                it = rlist->erase(it);
            else
                it++;
        }
    }

    // Add extra configurations not to be considered leading configurations.
    // These are added regardless of ConfigurationAverageEnergyRange
    int num_extra_configs = user_input.vector_variable_size("ExtraConfigurations");
    if(num_extra_configs)
    {
        pConfigList nrlist = std::make_shared<ConfigList>();
        nrlist->first.reserve(num_extra_configs);
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
            nrlist->first.push_back(extraconfig);
            nrlist->second++;
        }

        rlist->append(*GenerateRelativisticConfigurations(nrlist));
        rlist->SetSmallSize(rlist->size());
    }

    rlist->sort(FewestProjectionsFirstComparator());
    rlist->unique();
    rlist->SetSmallSize(rlist->size());

    // Print configuration list if requested
    if(user_input.search(2, "--print-relativistic-configurations", "--print-configurations"))
    {
        // Include energies
        if(one_body && two_body)
        {
            *outstream << "\nConfigurations; E; N:" << std::setprecision(6) << "\n";

            if(user_input.search("--print-relativistic-configurations"))
            {
                for(auto& config: *rlist)
                {
                    *outstream << config << "; "
                               << config.CalculateConfigurationAverageEnergy(orbitals->valence, one_body, two_body)
                               << "; " << config.GetNumberOfLevels() << "\n";
                }
            }
            else
            {
                pConfigList nrlist = GenerateNonRelConfigurations(rlist);
                for(auto& config: nrlist->first)
                {
                    *outstream << config << "; "
                               << config.CalculateConfigurationAverageEnergy(orbitals->valence, one_body, two_body)
                               << "; " << config.GetNumberOfLevels() << "\n";
                }
            }
            *outstream << std::flush;
        }
        else
        {
            *outstream << "\nConfigurations; N:" << std::endl;

            if(user_input.search("--print-relativistic-configurations"))
            {
                for(auto& config: *rlist)
                {
                    *outstream << config << "; " << config.GetNumberOfLevels() << "\n";
                }
            }
            else
            {
                pConfigList nrlist = GenerateNonRelConfigurations(rlist);
                for(auto& config: nrlist->first)
                {
                    *outstream << config << "; " << config.GetNumberOfLevels() << "\n";
                }
            }
            *outstream << std::flush;
        }
    }

    return rlist;
}

pRelativisticConfigList ConfigGenerator::GenerateRelativisticConfigurations(pConfigList nrlist) const
{
    pRelativisticConfigList rlist = std::make_shared<RelativisticConfigList>();
    auto itsmall = std::next(nrlist->first.begin(), nrlist->second);

    auto it = nrlist->first.begin();
    while(it != nrlist->first.end())
    {
        auto relconfs = it->relconfiglist;
        if(relconfs == nullptr)
            relconfs = it->GenerateRelativisticConfigs();

        rlist->append(*relconfs);
        if(it == itsmall)
            rlist->SetSmallSize(rlist->size());

        ++it;
    }

    rlist->sort(FewestProjectionsFirstComparator());
    rlist->unique();
    return rlist;
}

pConfigList ConfigGenerator::GenerateNonRelConfigurations(pRelativisticConfigListConst rlist) const
{
    pConfigList nrlist = std::make_shared<ConfigList>();
    nrlist->first.reserve(rlist->size());
    for(auto& rconf: *rlist)
    {
        nrlist->first.push_back(NonRelConfiguration(rconf));
    }
    SortAndUnique(nrlist);

    return nrlist;
}

pRelativisticConfigList ConfigGenerator::GenerateRelativisticConfigurations(pRelativisticConfigListConst rlist, const Symmetry& sym, pAngularDataLibrary angular_library) const
{
    pRelativisticConfigList prlist = std::make_shared<RelativisticConfigList>();
    RelativisticConfigList& returnedlist = *prlist;

    auto check_symmetry = [&sym](const RelativisticConfiguration& rconfig) {
        return (rconfig.GetParity() == sym.GetParity() && sym.GetTwoJ() <= rconfig.GetTwiceMaxProjection());
    };

    // [0, small_size)
    returnedlist = std::accumulate(rlist->begin(), rlist->small_end(), returnedlist,
                    [&check_symmetry](RelativisticConfigList& list, const RelativisticConfiguration& item) {
                        if(check_symmetry(item))
                            list.push_back(item);
                        return list;
                    });

    returnedlist.SetSmallSize(returnedlist.size());

    // [small_size, end)
    returnedlist = std::accumulate(rlist->small_end(), rlist->end(), returnedlist,
                    [&check_symmetry](RelativisticConfigList& list, const RelativisticConfiguration& item) {
                        if(check_symmetry(item))
                            list.push_back(item);
                        return list;
                    });

    if(angular_library)
        GenerateProjections(prlist, sym.GetTwoJ(), sym.GetTwoJ(), angular_library);

    return prlist;
}

void ConfigGenerator::GenerateExcitations(pRelativisticConfigList configlist, OrbitalMap& electron_valence, OrbitalMap& hole_valence) const
{
    // Go through the set of initial configurations
    RelativisticConfigList old_list(*configlist);

    auto config_it = old_list.begin();
    while(config_it != old_list.end())
    {
        // Move electrons and holes:
        // For each single particle state in the configuration
        auto particle_it = config_it->begin();
        while(particle_it != config_it->end())
        {
            // Electron or hole?
            if(particle_it->second > 0)
            {
                // Get another single particle state to move to
                for(const auto& electron: electron_valence)
                {
                    if(electron.first != particle_it->first)
                    {
                        RelativisticConfiguration new_config(*config_it);
                        new_config.RemoveSingleParticle(particle_it->first);
                        if(new_config.AddSingleParticle(electron.first))
                            configlist->push_back(new_config);
                    }
                }
            }
            else
            {   // Get another single particle state to move to
                for(const auto& hole: hole_valence)
                {
                    if(hole.first != particle_it->first)
                    {
                        RelativisticConfiguration new_config(*config_it);
                        new_config.AddSingleParticle(particle_it->first);
                        if(new_config.RemoveSingleParticle(hole.first))
                            configlist->push_back(new_config);
                    }
                }
            }

            particle_it++;
        }

        // Pair creation
        for(const auto& electron: electron_valence)
        {
            RelativisticConfiguration config_with_extra_electron(*config_it);
            if(config_with_extra_electron.AddSingleParticle(electron.first))
            {
                for(const auto& hole: hole_valence)
                {
                    RelativisticConfiguration new_config(config_with_extra_electron);
                    if(new_config.RemoveSingleParticle(hole.first))
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
                    RelativisticConfiguration new_config(*config_it);
                    new_config.AddSingleParticle(particle_it->first);
                    new_config.RemoveSingleParticle(other_particle_it->first);
                    configlist->push_back(new_config);
                }
                else if((particle_it->second > 0) && (other_particle_it->second < 0))
                {
                    RelativisticConfiguration new_config(*config_it);
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

    configlist->SetSmallSize(configlist->size());
    configlist->sort();
    configlist->unique();
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
