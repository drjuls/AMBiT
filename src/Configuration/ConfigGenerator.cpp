#include "Include.h"
#include "ConfigGenerator.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Specification/Specification.h"

namespace Ambit
{
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

pRelativisticConfigList ConfigGenerator::GenerateConfigurations(pHFIntegrals one_body, pSlaterIntegrals two_body)
{
    auto config = GlobalSpecification::Instance();
    // Do sublist first since leading_configs should come from CI, rather than CI/SmallSide
    pRelativisticConfigList sublist = nullptr;
    if(config->ci_do_smallside)
    {
        sublist = ParseAndGenerateConfigurations(&(config->ci_smallside), one_body, two_body);
    }

    pRelativisticConfigList biglist = ParseAndGenerateConfigurations(&(config->ci),one_body, two_body);

    // Join sublist and biglist by getting set difference and appending to sublist.
    if(sublist)
    {
        pRelativisticConfigList setdiff = std::make_shared<RelativisticConfigList>(biglist->size());
        auto last = std::set_difference(biglist->begin(), biglist->end(), sublist->begin(), sublist->end(), setdiff->begin());
        setdiff->erase(last, setdiff->end());

        biglist = sublist;
        biglist->append(*setdiff);
    }

    return biglist;
}

pRelativisticConfigList ConfigGenerator::ParseAndGenerateConfigurations(CIBaseSpecification* ci_spec, pHFIntegrals one_body, pSlaterIntegrals two_body)
{
    pRelativisticConfigList rlist = std::make_shared<RelativisticConfigList>();

    int electron_excitations = ci_spec->electron_excitations;
    int hole_excitations = ci_spec->hole_excitations;
    
    // Get leading configurations
    leading_configs->first.clear();
    int numValenceElectrons = 0;
    
    // We must have at least one leading configuration for a well-formed calculation. All
    // other "leading config" or "extra config" inputs are optional
    unsigned num_configs = ci_spec->leading_configurations.size();
    auto spec_lcs = ci_spec->leading_configurations;
    leading_configs->first.reserve(num_configs);
    for(unsigned i = 0; i < num_configs; i++)
    {
        const std::string name = spec_lcs[i];
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

    leading_configs->second = leading_configs->first.size();
    SortAndUnique(leading_configs);
    rlist->append(*GenerateRelativisticConfigurations(leading_configs));

    // Get any leading relativistic configurations

    if(ci_spec->leading_rel_configurations)
    {
        auto spec_rlcs = ci_spec->leading_rel_configurations.value();
        num_configs = spec_rlcs.size();

        for(unsigned i = 0; i < num_configs; i++)
        {
            const std::string name = spec_rlcs[i];
            RelativisticConfiguration config(name);

            // Check that the configuration gels with the number of electrons
            if(leading_configs->first.size() == 0)
                numValenceElectrons = config.ElectronNumber();
            else if(config.ElectronNumber() != numValenceElectrons)
            {   *errstream << "USAGE: LeadingConfiguration " << name
                           << " does not have correct number of valence electrons." << std::endl;
                exit(1);
            }
            rlist->push_back(config);
            leading_configs->first.push_back(NonRelConfiguration(config));
        }

        leading_configs->second = leading_configs->first.size();
        SortAndUnique(leading_configs);
    }

    // Get maxenergy
    bool use_maxenergy = false;
    double upper_energy = 0.0;
    if(one_body && two_body && ci_spec->configuration_average_energy_range)
    {
        use_maxenergy = true;
        upper_energy = ci_spec->configuration_average_energy_range.value().second;
    }

    // Total number of excitation steps
    int total_num_excitations = mmax(electron_excitations, hole_excitations);

    for(int excitation_step = 0; excitation_step < total_num_excitations; excitation_step++)
    {
        OrbitalMap valence_electrons(orbitals->GetLattice());
        OrbitalMap valence_holes(orbitals->GetLattice());

        // Electron subset
        if(excitation_step < electron_excitations)
        {
            valence_electrons = *orbitals->particle;

            if(ci_spec->electron_bounds)
            {
                auto bounds = ci_spec->electron_bounds.value();
                std::string CI_basis_string = bounds[excitation_step];
                // Limits here are in the form of ValenceBasis (pqn up to and including limit)
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
            valence_holes = *orbitals->hole;

            if(ci_spec->hole_bounds)
            {
                auto bounds = ci_spec->hole_bounds.value();
                std::string CI_basis_string = bounds[excitation_step];
                // Limits here are in the form of FrozenCore (pqn should be bigger than limit)
                std::vector<int> limits = ConfigurationParser::ParseBasisSize(CI_basis_string);
                
                auto valence_it = valence_holes.begin();
                while(valence_it != valence_holes.end())
                {
                    int L = valence_it->first.L();
                    
                    if(L < limits.size() && valence_it->first.PQN() <= limits[L])
                        valence_it = valence_holes.erase(valence_it);
                    else
                        valence_it++;
                }
            }
        }

        GenerateExcitations(rlist, valence_electrons, valence_holes);

        if(use_maxenergy && excitation_step < total_num_excitations-1)
        {   auto newend = std::remove_if(rlist->begin(), rlist->end(),
                            [&](RelativisticConfiguration& item){
                                return (item.CalculateConfigurationAverageEnergy(orbitals->valence, one_body, two_body) > upper_energy);
                            });
            rlist->erase(newend, rlist->end());
        }
    }

    // Trim configurations outside of ConfigurationAverageEnergyRange
    if(one_body && two_body && ci_spec->configuration_average_energy_range)
    {
        double lower_energy = ci_spec->configuration_average_energy_range.value().first;

        auto newend = std::remove_if(rlist->begin(), rlist->end(),
                        [&](RelativisticConfiguration& item){
                            double energy = item.CalculateConfigurationAverageEnergy(orbitals->valence, one_body, two_body);
                            return (energy < lower_energy || energy > upper_energy);
                        });
        rlist->erase(newend, rlist->end());
    }

    // Add extra configurations not to be considered leading configurations.
    // These are added regardless of ConfigurationAverageEnergyRange
    if(ci_spec->extra_configurations)
    {
        auto spec_extra_configs = ci_spec->extra_configurations.value();
        unsigned num_extra_configs = spec_extra_configs.size();

        pConfigList nrlist = std::make_shared<ConfigList>();
        nrlist->first.reserve(num_extra_configs);
        for(unsigned i = 0; i < num_extra_configs; i++)
        {
            const std::string extraname = spec_extra_configs[i];
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

        // Now add to leading configurations for sigma 3
        leading_configs->first.insert(leading_configs->first.end(), nrlist->first.begin(), nrlist->first.end());
    }

    // Add extra relativistic configurations not to be considered leading configurations for the purpose of configuration generation.
    if(ci_spec->extra_rel_configurations)
    {

        auto spec_extra_rconfigs = ci_spec->extra_rel_configurations.value();
        unsigned num_extra_configs = spec_extra_rconfigs.size();

        for(unsigned i = 0; i < num_extra_configs; i++)
        {
            const std::string extraname = spec_extra_rconfigs[i];
            RelativisticConfiguration extraconfig(extraname);

            if(extraconfig.ElectronNumber() != numValenceElectrons)
            {
                *errstream << "USAGE: ExtraConfiguration " << extraname
                           << " does not have correct number of valence electrons." << std::endl;
                exit(1);
            }
            rlist->push_back(extraconfig);
        }
    }

    rlist->SetSmallSize(rlist->size());
    rlist->sort(MostWorkFirstComparator());
    rlist->unique();
    rlist->SetSmallSize(rlist->size());

    // Print configuration list if requested
    if(ci_spec->print_configurations || ci_spec->print_rel_configurations)
    {
        // Include energies
        if(one_body && two_body)
        {
            *outstream << "\nConfigurations; E; N:" << std::setprecision(6) << "\n";

            if(ci_spec->print_rel_configurations)
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
            *outstream << std::endl;
        }
        else
        {
            *outstream << "\nConfigurations; N:" << std::endl;

            if(ci_spec->print_rel_configurations)
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
            *outstream << std::endl;
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
    if(it == itsmall)
        rlist->SetSmallSize(rlist->size());

    rlist->sort(MostWorkFirstComparator());
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
    nrlist->second = rlist->small_size();
    SortAndUnique(nrlist);

    return nrlist;
}

pRelativisticConfigList ConfigGenerator::GenerateRelativisticConfigurations(pRelativisticConfigListConst rlist, const Symmetry& sym, pAngularDataLibrary angular_library) const
{
    auto config = GlobalSpecification::Instance();
    pRelativisticConfigList prlist = std::make_shared<RelativisticConfigList>();
    RelativisticConfigList& returnedlist = *prlist;

    auto check_symmetry = [&sym](const RelativisticConfiguration& rconfig) {
        return (rconfig.GetParity() == sym.GetParity() && sym.GetTwoJ() <= rconfig.GetTwiceMaxProjection());
    };

    /* Select the configurations from the main list with the desired parity
     * NOTE: we write this loop out by hand rather than using std::accumulate since accumulate casts its
     * arguments to const, which causes lots of copy allocations (even when passing arguments by 
     * reference). This drastically reduces performance when rlist is large. 
     */
    //[0, small_size)
    auto it = rlist->begin();
    while(it != rlist->small_end())
    {
        if(check_symmetry(*it))
        {
            returnedlist.push_back(*it);
        }
        it++;
    }

    returnedlist.SetSmallSize(returnedlist.size());

    // [small_size, end)
    it = rlist->small_end();
    while(it != rlist->end())
    {
        if(check_symmetry(*it))
        {
            returnedlist.push_back(*it);
        }
        it++;
    }

    if(angular_library)
        GenerateProjections(prlist, sym, sym.GetTwoJ(), angular_library);

    //if(user_input.search("CI/--sort-matrix-by-configuration"))
    if(config->ci.sort_matrix_by_configuration)
        returnedlist.sort(ConfigurationComparator());
    else
        returnedlist.sort(MostWorkFirstComparator());

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
    configlist->sort(MostProjectionsFirstComparator());
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
}
