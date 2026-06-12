#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif

#ifdef AMBIT_USE_OPENMP
#include<omp.h>
#endif

#include "Include.h"
#include "Atom.h"
#include "Universal/Enums.h"
#include "HartreeFock/NonRelInfo.h"
#include "Configuration/ConfigGenerator.h"
#include "Configuration/HamiltonianMatrix.h"
#include "Configuration/GFactor.h"
#include "MBPT/OneElectronMBPT.h"
#include "MBPT/CoreValenceIntegrals.h"
#include "MBPT/BruecknerDecorator.h"
#include "HartreeFock/HartreeFocker.h"
#include "HamiltonianTypes.h"
#include "Universal/Profiler.h"

namespace Ambit
{
void Atom::MakeMBPTIntegrals()
{
    auto specification = GlobalSpecification::Instance();
    bool check_sizes = specification->check_sizes;
    bool one_body_mbpt = specification->mbpt_one_body;
    bool two_body_mbpt = specification->mbpt_two_body;

    // Bare integrals for MBPT
    pHFIntegrals bare_one_body_integrals = std::make_shared<HFIntegrals>(orbitals, hf);
    pSlaterIntegrals bare_two_body_integrals = std::make_shared<SlaterIntegralsFlatHash>(orbitals, hartreeY);

    // MBPT calculators
    std::optional<std::string> fermi_orbitals = specification->mbpt_energy_denom_orbitals;
    pCoreMBPTCalculator core_mbpt = std::make_shared<CoreMBPTCalculator>(orbitals, bare_one_body_integrals, bare_two_body_integrals, fermi_orbitals);
    pValenceMBPTCalculator val_mbpt = std::make_shared<ValenceMBPTCalculator>(orbitals, bare_one_body_integrals, bare_two_body_integrals, fermi_orbitals);

    auto& valence = orbitals->valence;

    pOneElectronMBPT mbpt_integrals_one = std::make_shared<OneElectronMBPT>(orbitals, hf, core_mbpt, val_mbpt, identifier + ".one.int");
    pCoreValenceIntegralsMap mbpt_integrals_two = std::make_shared<CoreValenceIntegralsMap>(orbitals, core_mbpt, val_mbpt, identifier + ".two.int");

    // Use subtraction diagrams and extra box diagrams?
    // First find out whether the core is our closed shell core, while we're at it, get maximum pqn.
    int max_pqn_in_core = 0;
    bool is_open_shell = (orbitals->core->size() != open_core->size()); // Minimum requirement
    for(auto pair: *orbitals->core)
    {
        double occ = open_core->GetOccupancy(pair.first);
        if(fabs(occ - double(pair.first.MaxNumElectrons())) > 0.01)
            is_open_shell = true;

        max_pqn_in_core = mmax(max_pqn_in_core, pair.first.PQN());
    }

    // We can always force subtraction diagrams, and it's always a good idea if there have been injected orbitals
    if(specification->mbpt_force_subtraction || specification->basis_inject_orbitals)
        is_open_shell = true;
    if(specification->mbpt_no_subtraction)
        is_open_shell = false;

    bool use_box = !specification->mbpt_no_extra_box;
    bool include_core = !specification->mbpt_no_core;
    mbpt_integrals_one->IncludeCore(include_core, include_core && is_open_shell);
    mbpt_integrals_two->IncludeCore(include_core, include_core && is_open_shell, include_core && use_box);

    bool include_valence = specification->mbpt_use_valence;
    mbpt_integrals_one->IncludeValence(include_valence && is_open_shell);
    mbpt_integrals_two->IncludeValence(include_valence, include_valence && is_open_shell, include_valence && use_box);

    // Adjust delta
    double delta = specification->mbpt_delta;
    core_mbpt->SetEnergyShift(delta);
    val_mbpt->SetEnergyShift(delta);

    // Also set a floor for the valence-MBPT energy denominators
    double energy_denom_floor = specification->mbpt_energy_denom_floor;
    *logstream << "MBPT denominator floor = " << energy_denom_floor << std::endl;
    core_mbpt->SetEnergyFloor(energy_denom_floor);
    val_mbpt->SetEnergyFloor(energy_denom_floor);

    // Calculate two electron integrals on valence orbitals with limits on the PQNs of the orbitals.
    // Use max_pqn_1, max_pqn_2 and max_pqn_3 to keep size down.
    // For two electron integrals:
    //      i1.pqn <= limit1
    //      (i2.pqn or i3.pqn) <= limit2
    //      (i2.pqn and i3.pqn) <= limit3
    // For 'x'* 3 waves (spd) and 'y'* 4 waves (spdf) in basis set,
    //      N = 61 x^4       N = 279 y^4.
    // After max_pqn_1 = 4,
    //      N ~ 502 x^3      N ~ 1858 y^3,
    // and hopefully after max_pqn_2 and then max_pqn_3
    //      N~ x^2 and then N ~ x, respectively.
    std::vector<pOrbitalMap> valence_subset(4, valence);

    if(two_body_mbpt || check_sizes)
    {
        if(specification->mbpt_twobody_storage_limits)
        {
            std::vector<unsigned int> limits = specification->mbpt_twobody_storage_limits.value();
            unsigned int num_limits = limits.size();
            for(unsigned int i = 0; i < mmin(num_limits, 4); i++)
            {
                int max_pqn = limits[i];

                if(max_pqn)
                {   // Change orbital
                    valence_subset[i].reset(new OrbitalMap(*valence));
                    auto it = valence_subset[i]->begin();
                    while(it != valence_subset[i]->end())
                    {
                        if(it->first.PQN() > max_pqn)
                            it = valence_subset[i]->erase(it);
                        else
                            it++;
                    }
                }
            }
        }
        else
        {   // Default when not specified
            int max_pqn = max_pqn_in_core + 2;
            valence_subset[0].reset(new OrbitalMap(*valence));
            auto it = valence_subset[0]->begin();
            while(it != valence_subset[0]->end())
            {
                if(it->first.PQN() > max_pqn)
                    it = valence_subset[0]->erase(it);
                else
                    it++;
            }

            valence_subset[1] = valence_subset[0];
        }
    }

    if(check_sizes)
    {
        unsigned int size = 0;
        if(include_core)
            size += core_mbpt->GetStorageSize();
        if(include_valence)
            size += val_mbpt->GetStorageSize();

        *outstream << "Num stored coulomb integrals for MBPT = " << size;

        size = mbpt_integrals_one->CalculateOneElectronIntegrals(valence, valence, true);
        *outstream << "\nNum one-body mbpt integrals: " << size;

        size = mbpt_integrals_two->CalculateTwoElectronIntegrals(valence_subset[0], valence_subset[1], valence_subset[2], valence_subset[3], true);
        *outstream << "\nNum two-body mbpt integrals: " << size << std::endl;
    }
    else
    {
        if(one_body_mbpt)
        {
            // Start the timer
            auto profiler = Profiler::Instance();
            profiler->one_body_mbpt.start();

            mbpt_integrals_one->Read(identifier + ".one.int");
            unsigned int size = mbpt_integrals_one->CalculateOneElectronIntegrals(valence, valence);

            // Stop the timer
            profiler->one_body_mbpt.stop();
            *logstream << "One-body MBPT integrals complete: size = " << size << std::endl;
        }

        if(two_body_mbpt)
        {
            // Start the timer
            auto profiler = Profiler::Instance();
            profiler->two_body_mbpt.start();

            mbpt_integrals_two->Read(identifier + ".two.int");
            unsigned int size = mbpt_integrals_two->CalculateTwoElectronIntegrals(valence_subset[0], valence_subset[1], valence_subset[2], valence_subset[3]);

            // Stop the timer
            profiler->two_body_mbpt.stop();
            *logstream << "Two-body MBPT integrals complete: size = " << size << std::endl;
        }
    }
}

void Atom::MakeCIIntegrals()
{
    // Start the timer
    auto profiler = Profiler::Instance();
    profiler->slater.start();

    auto specification = GlobalSpecification::Instance();
    ClearIntegrals();

    pSlaterIntegrals two_body_integrals;
    bool one_body_mbpt = specification->mbpt_one_body;
    bool two_body_mbpt = specification->mbpt_two_body;
    bool three_body_mbpt = specification->mbpt_three_body;

    if(!two_body_mbpt)
        two_body_integrals.reset(new SlaterIntegralsFlatHash(orbitals, hartreeY));
    else
        two_body_integrals.reset(new SlaterIntegralsFlatHash(orbitals, hartreeY, false));

    if(three_body_mbpt)
    {
        std::optional<std::string> fermi_orbitals = specification->mbpt_energy_denom_orbitals;
        threebody_electron = std::make_shared<Sigma3Calculator>(orbitals, two_body_integrals, fermi_orbitals);
        threebody_electron->IncludeCore(!specification->mbpt_no_core);
        threebody_electron->IncludeValence(specification->mbpt_use_valence);
    }

    auto& valence = orbitals->valence;

    // We don't need the two body integrals if we are just checking sizes and they are not needed by ConfigGenerator.
    if(specification->check_sizes
       && !specification->ci.configuration_average_energy_range
       && !specification->ci_smallside.configuration_average_energy_range
       && !specification->ci.print_configurations && !specification->ci.print_rel_configurations
       && !specification->ci_smallside.print_configurations && !specification->ci_smallside.print_rel_configurations)
    {
        unsigned int size = two_body_integrals->CalculateTwoElectronIntegrals(valence, valence, valence, valence, true);
        *outstream << "\nNum Coulomb integrals: " << size << std::endl;

        if(three_body_mbpt)
            *outstream << "\nSigma3 Coulomb integrals: " << threebody_electron->GetStorageSize() << std::endl;
    }
    else
    {   // Normal integrals
        hf_electron.reset(new HFIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(valence, valence);

        // Don't need two body integrals if we're not doing CI
        if(!specification->no_ci)
        {
            // Check if we can read existing Slater Integrals from disk
            bool use_read = !specification->clean_run;

            // Generate filename
            std::string slater_filename = identifier + ".coulomb.int";
            // Convert spaces to underscores in filename
            std::replace_if(slater_filename.begin(), slater_filename.end(),
                            [](char c){ return (c =='\r' || c =='\t' || c == ' ' || c == '\n');}, '_');

            // Read Slater integrals if available, but only if we don't have the "-c" flag
            bool generate_integrals = true;

            // If reading from a checkpoint file, make sure the file exists and is
            // valid, then read
            if(use_read)
            {
                generate_integrals = !two_body_integrals->Read(slater_filename);
                if(!generate_integrals)
                    *logstream << "Read Slater integrals from file "
                               << slater_filename << std::endl;
            }

            // Generate the integrals if reading from file failed, or if we're doing a
            // clean run
            if(generate_integrals)
            {
                *logstream << "Generating Slater integrals..." << std::endl;
                two_body_integrals->CalculateTwoElectronIntegrals(valence, valence, valence, valence);

                // And write the resulting integrals to the file
                two_body_integrals->Write(slater_filename);
            }

            // Print size of integrals if requested
            if(specification->check_sizes)
                *outstream << "\nNum Coulomb integrals: " << two_body_integrals->size() << std::endl;
        }

        // Add stored MBPT integrals
        if(one_body_mbpt)
        {
            if(specification->mbpt_onebody_scaling)
            {
                std::vector<double> scaling = specification->mbpt_onebody_scaling.value();
                std::map<int, double> scaling_map;
                for(int i = 0; i < scaling.size()-1; i+=2)
                {
                    int kappa = scaling[i];
                    double lambda = scaling[i+1];

                    scaling_map[kappa] = lambda;
                }

                hf_electron->Read(identifier + ".one.int", scaling_map);
            }
            else
                hf_electron->Read(identifier + ".one.int");
        }
        if(two_body_mbpt)
            two_body_integrals->Read(identifier + ".two.int");

        bool include_box = two_body_mbpt && !specification->mbpt_no_extra_box;
        bool include_off_parity = include_box || hartreeY->OffParityExists();
        twobody_electron = std::make_shared<TwoElectronCoulombOperator>(two_body_integrals, include_off_parity);

        // Calculate integrals for sigma3
        if(three_body_mbpt)
        {
            threebody_electron->UpdateIntegrals();
            if(specification->check_sizes)
                *outstream << "\nSigma3 Coulomb integrals: " << threebody_electron->GetStorageSize() << std::endl;
        }
    }

    // Stop the timer
    profiler->slater.stop();
}

void Atom::ClearIntegrals()
{
    hf_electron = nullptr;
    twobody_electron = nullptr;
    threebody_electron = nullptr;
}

void Atom::InitialiseAngularDataLibrary(pAngularDataLibrary trial)
{
    if(trial)
        angular_library = trial;
    else if(angular_library == nullptr)
    {
        std::string angular_directory = string_macro(ANGULAR_DATA_DIRECTORY);
        angular_library = std::make_shared<AngularDataLibrary>(angular_directory);
    }
}

pLevelStore Atom::ChooseHamiltoniansAndRead(pAngularDataLibrary angular_lib)
{
    auto specification = GlobalSpecification::Instance();
    // Read existing levels?
    bool use_read = !specification->clean_run;

    if(specification->configuration_average)
    {
        ConfigGenerator gen(orbitals);
        if(twobody_electron == nullptr)
            MakeCIIntegrals();
        allconfigs = gen.GenerateConfigurations(hf_electron, twobody_electron->GetIntegrals());

        // Create empty level set
        levels = std::make_shared<LevelMap>(identifier, nullptr);
        return levels;
    }

    // Get angular library
    InitialiseAngularDataLibrary(angular_lib);

    if(specification->ci.memory_saver)
    {
        // LevelDirectory must have a value if we're running with --memory-saver and this has been
        // enforced by the specification
        std::string dir = specification->level_directory.value();
        if(dir.empty())
            levels = std::make_shared<FileSystemLevelStore>(identifier, angular_library);
        else
            levels = std::make_shared<FileSystemLevelStore>(dir, identifier, angular_library);
    }
    else if(use_read)
    {
        pHamiltonianID key;
        if(specification->no_ci)
            key = std::make_shared<SingleOrbitalID>();
        else if(specification->ci.single_configuration_ci)
            key = std::make_shared<NonRelID>();
        else
            key = std::make_shared<HamiltonianID>();

        levels = std::make_shared<LevelMap>(key, identifier, angular_library);
    }
    else
        levels = std::make_shared<LevelMap>(identifier, angular_library);

    if(specification->no_ci)
    {
        // Use all symmetries from valence set
        pHamiltonianID key;
        for(auto& it: *orbitals->valence)
        {
            key = std::make_shared<SingleOrbitalID>(it.first);
            levels->keys.insert(key);
        }
    }
    else
    {   // Generate non-rel configurations and hence choose Hamiltonians
        ConfigGenerator gen(orbitals);

        if(twobody_electron == nullptr)
            MakeCIIntegrals();

        if(twobody_electron)
            allconfigs = gen.GenerateConfigurations(hf_electron, twobody_electron->GetIntegrals());
        else
            allconfigs = gen.GenerateConfigurations();

        leading_configs = gen.GetLeadingConfigs();

        ChooseHamiltonians(allconfigs);
    }

    return levels;
}

pLevelStore Atom::ChooseHamiltonians(pRelativisticConfigList rlist)
{
    auto specification = GlobalSpecification::Instance();
    // TODO EVK: the user input has these as unsigned integers (they have to be, because it doesn't
    // make any sense to pass a negative value for 2J in the config file), but all the
    // HamiltonianID objects take signed integers. Should probably change this over, but would need
    // to modify the binary file formats for checkpoint files.
    std::vector<unsigned int> even_symmetries;
    std::vector<unsigned int> odd_symmetries;
    pHamiltonianID key;

    // Get user symmetries if CI/--all-symmetries is not used. Fine to just do a deep-copy from the
    // "specification" struct
    if(!specification->ci.all_symmetries)
    {
        // Even parity
        even_symmetries = specification->ci.even_parity_twoj;
        // Odd parity
        odd_symmetries = specification->ci.odd_parity_twoj;
    }

    // Single configuration CI: each Hamiltonian only uses one non-rel configuration
    if(specification->ci.single_configuration_ci)
    {
        ConfigGenerator gen(orbitals);
        pConfigList nrlist = gen.GenerateNonRelConfigurations(rlist);

        if(specification->ci.all_symmetries)
        {
            // Populate levels with all symmetries
            for(auto& nrconfig: nrlist->first)
            {
                int maxTwoJ = nrconfig.GetTwiceMaxProjection();
                for(int two_j = maxTwoJ%2; two_j <= maxTwoJ; two_j += 2)
                {
                    key = std::make_shared<NonRelID>(nrconfig, two_j);
                    levels->keys.insert(key);
                }
            }
        }
        else
        {   // Populate levels with symmetries found in sets.
            for(auto& nrconfig: nrlist->first)
            {
                Parity P = nrconfig.GetParity();
                unsigned int maxTwoJ = nrconfig.GetTwiceMaxProjection();

                if(P == Parity::even)
                {
                    for(unsigned int& two_j: even_symmetries)
                        if(two_j <= maxTwoJ)
                        {
                            key = std::make_shared<NonRelID>(nrconfig, two_j);
                            levels->keys.insert(key);
                        }
                }
                else
                {   for(unsigned int& two_j: odd_symmetries)
                        if(two_j <= maxTwoJ)
                        {
                            key = std::make_shared<NonRelID>(nrconfig, two_j);
                            levels->keys.insert(key);
                        }
                }
            }
        }
    }
    // Standard CI: all configs together
    else
    {   // Get symmetries
        if(specification->ci.all_symmetries)
        {
            int max_twoJ_even = -1;
            int max_twoJ_odd = -1;

            for(auto& config: *rlist)
            {
                if(config.GetParity() == Parity::even)
                    max_twoJ_even = mmax(max_twoJ_even, config.GetTwiceMaxProjection());
                else
                    max_twoJ_odd = mmax(max_twoJ_odd, config.GetTwiceMaxProjection());
            }

            for(unsigned int twoJ = mmax(max_twoJ_even%2, 0); twoJ <= max_twoJ_even; twoJ+=2)
                even_symmetries.push_back(twoJ);
            for(unsigned int twoJ = mmax(max_twoJ_odd%2, 0); twoJ <= max_twoJ_odd; twoJ+=2)
                odd_symmetries.push_back(twoJ);
        }

        // Populate levels
        for(unsigned int& two_j: even_symmetries)
        {   key = std::make_shared<HamiltonianID>(two_j, Parity::even);
            levels->keys.insert(key);
        }
        for(unsigned int& two_j: odd_symmetries)
        {   key = std::make_shared<HamiltonianID>(two_j, Parity::odd);
            levels->keys.insert(key);
        }
    }

    if(levels->keys.empty())
    {   *errstream << "USAGE: No symmetries requested (EvenParityTwoJ or OddParityTwoJ)" << std::endl;
    }

    return levels;
}

/** Check sizes of matrices before doing full scale calculation. */
void Atom::CheckMatrixSizes(pAngularDataLibrary angular_lib)
{
    auto specification = GlobalSpecification::Instance();
    if(specification->no_ci)
        return;

    // Get angular library
    InitialiseAngularDataLibrary(angular_lib);

    // Don't read existing levels
    levels = std::make_shared<LevelMap>(identifier, angular_library);

    // CI integrals
    MakeCIIntegrals();

    // Generate configurations again; don't read from disk. */
    ConfigGenerator gen(orbitals);
    if(twobody_electron)
        allconfigs = gen.GenerateConfigurations(hf_electron, twobody_electron->GetIntegrals());
    else
        allconfigs = gen.GenerateConfigurations();
    leading_configs = gen.GetLeadingConfigs();

    ChooseHamiltonians(allconfigs);

    // Get complete list of relativistic configs for all symmetries
    // to help with load balancing for AngularData to create CSFs.
    std::vector<Symmetry> symmetries;
    for(auto& key: levels->keys)
    {
        symmetries.emplace_back(key->GetSymmetry());
    }
    std::sort(symmetries.begin(), symmetries.end());
    auto end = std::unique(symmetries.begin(), symmetries.end());
    symmetries.erase(end, symmetries.end());

    // Generate CSFs
    *outstream << "Calculating CSFs: " << std::endl;
    unsigned int total_levels = 0;

    // Do all configs in a symmetry at once
    for(auto& sym: symmetries)
    {
        // Get list of relativistic configs
        pRelativisticConfigList configs = gen.GenerateRelativisticConfigurations(allconfigs, sym);

        *outstream << "J(P) = " << sym.GetJ() << "(" << ShortName(sym.GetParity()) << "): "
                   << std::setw(6) << std::right << configs->size();
        if(configs->small_size() < configs->size())
            *outstream << " x " << std::setw(6) << configs->small_size();
        *outstream << " rel. configurations; " << std::flush;

        // Generate all projections for this symmetry and write
        gen.GenerateProjections(configs, sym, sym.GetTwoJ(), angular_library);

        *outstream << configs->NumCSFs();
        if(configs->small_size() < configs->size())
            *outstream << " x " << configs->NumCSFsSmall();
        *outstream << " CSFs." << std::endl;

        total_levels += configs->NumCSFs();
    }

    *outstream << "\nTotal number of levels (all symmetries included) = " << total_levels << std::endl;

    // Get Hamiltonian sizes for single configurations
    if(specification->ci.single_configuration_ci)
    {
        *outstream << "\nHamiltonian matrix sizes: " << std::endl;
        for(auto& key: levels->keys)
        {
            pConfigList nrconfiglist = std::make_shared<ConfigList>();
            nrconfiglist->first.emplace_back(std::static_pointer_cast<const NonRelID>(key)->GetNonRelConfiguration());
            nrconfiglist->second = 1;
            auto rconfigs = gen.GenerateRelativisticConfigurations(nrconfiglist);

            unsigned int num_CSFs = 0;
            for(auto& config: *rconfigs)
            {
                num_CSFs += angular_library->GetData(config, key->GetSymmetry(), key->GetTwoJ())->NumCSFs();
            }

            *outstream << key->Print() << ":"
                       << std::setw(6) << std::right << rconfigs->size() << " rel. configurations; "
                       << num_CSFs <<  " CSFs." << std::endl;
        }
    }
}

pLevelStore Atom::CalculateEnergies()
{
    ChooseHamiltoniansAndRead();

    for(auto& key: *levels)
        CalculateEnergies(key);

    return levels;
}

LevelVector Atom::CalculateEnergies(pHamiltonianID hID)
{
    auto specification = GlobalSpecification::Instance();
    // Check whether to reuse levels and matrix files
    bool use_read = !specification->clean_run;

    // This function is public and can call the other CalculateEnergies variants.
    LevelVector levelvec = levels->GetLevels(hID);
    std::string filename = identifier + ".levels";

    auto nrID = std::dynamic_pointer_cast<NonRelID>(hID);
    auto soID = std::dynamic_pointer_cast<SingleOrbitalID>(hID);

    if(soID)
    {
        if(levelvec.NumLevels() == 0)
            levelvec = SingleElectronConfigurations(hID);
    }
    else
    {   // Get relativistic configurations
        pRelativisticConfigList& configs = levelvec.configs;

        if(configs == nullptr)
        {
            ConfigGenerator gen(orbitals);

            if(nrID)
            {
                pConfigList nrconfiglist = std::make_shared<ConfigList>();
                nrconfiglist->first.emplace_back(nrID->GetNonRelConfiguration());
                nrconfiglist->second = 1;
                configs = gen.GenerateRelativisticConfigurations(nrconfiglist);
                configs = gen.GenerateRelativisticConfigurations(configs, nrID->GetSymmetry(), angular_library);
            }
            else
            {
                if(configs == nullptr)
                {
                    if(allconfigs == nullptr)
                    {
                        if(specification->ci.configuration_average_energy_range
                           || specification->ci_smallside.configuration_average_energy_range
                           || specification->ci.print_configurations && !specification->ci.print_rel_configurations
                           || specification->ci_smallside.print_configurations && !specification->ci_smallside.print_rel_configurations)
                        {
                            if(twobody_electron == nullptr)
                                MakeCIIntegrals();

                            allconfigs = gen.GenerateConfigurations(hf_electron, twobody_electron->GetIntegrals());
                        }
                        else
                            allconfigs = gen.GenerateConfigurations();

                        leading_configs = gen.GetLeadingConfigs();
                    }
                }
                configs = gen.GenerateRelativisticConfigurations(allconfigs, hID->GetSymmetry(), angular_library);
            }

            angular_library->RemoveUnused();
        }

        // Only continue if we don't have enough levels
        int num_solutions = specification->ci.num_solutions;
        num_solutions = (num_solutions? mmin(num_solutions, configs->NumCSFs()): configs->NumCSFs());
        bool do_CI = (levelvec.NumLevels() < num_solutions);
        if(do_CI)
        {
            if(twobody_electron == nullptr)
                MakeCIIntegrals();

            // If we're using OpenMP then the chunksize should be a multiple of the number of threads
            int chunksize = specification->ci.chunksize;

            std::unique_ptr<HamiltonianMatrix> H;
            if(threebody_electron)
                H.reset(new HamiltonianMatrix(hf_electron, twobody_electron, threebody_electron, leading_configs, configs));
            else
                H.reset(new HamiltonianMatrix(hf_electron, twobody_electron, configs));

            // Generate filename
            std::string hamiltonian_filename = identifier + "." + hID->Name() + ".matrix";
            // Convert spaces to underscores in filename
            std::replace_if(hamiltonian_filename.begin(), hamiltonian_filename.end(),
                            [](char c){ return (c =='\r' || c =='\t' || c == ' ' || c == '\n');}, '_');

            // Read Hamiltonian if available, but only if we don't have the "-c" flag
            bool generate_matrix = true;
            // Do not solve the matrix, just generate Hamiltonian and write to file
            bool dont_solve_matrix = specification->no_diagonalisation;

            // If reading from a checkpoint file, make sure the file exists and is valid, then
            // read
            if(use_read)
            {
                generate_matrix = !H->Read(hamiltonian_filename, chunksize);
            }

            // Only generate if we're doing a clean run or if reading from the matrix checkpoint
            // file failed for some reason
            if(generate_matrix)
            {
                H->GenerateMatrix(chunksize);
                //H->PollMatrix();

                if(specification->ci.output_write_hamiltonian || dont_solve_matrix)
                    H->Write(hamiltonian_filename);
            }

            if(specification->ci.output_print_hamiltonian)
            {
                auto rel_it = configs->begin();
                while(rel_it != configs->end())
                {
                    *outstream << rel_it->Name();
                    if(rel_it++ != configs->end())
                    {
                        *outstream << ",";
                    }
                }
                *outstream << std::endl;

                *outstream << std::setprecision(12);
                *outstream << "Matrix Before:\n" << *H << std::endl;
            }

            if(dont_solve_matrix)
            {   *outstream << "\n";
                return {};  // Return empty levelvector
            }

            #ifdef AMBIT_USE_SCALAPACK
            if(specification->ci.scalapack)
            {
                // Write matrix if it hasn't already been written
                if(!user_input.search("CI/Output/--write-hamiltonian"))
                    H->Write(hamiltonian_filename);

                if(specification->ci.max_energy)
                {
                    double max_energy = specification->ci.max_energy.value());
                    levelvec = H->SolveMatrixScalapack(hID, configs->NumCSFs(), hamiltonian_filename, max_energy);
                }
                else
                {
                    levelvec = H->SolveMatrixScalapack(hID, num_solutions, hamiltonian_filename);
                }
            }
            else
            #endif
            levelvec = H->SolveMatrix(hID, num_solutions, hamiltonian_filename);
            levels->Store(hID, levelvec);
        }

        // Check if gfactor overrides are present, otherwise decide on course of action
        bool get_gfactors = levels->GFactorsNeeded();

        // Don't bother doing any checks if we've got pre-calculated g-factors stored
        if(get_gfactors)
        {
            if(hID->GetTwoJ() == 0 || specification->ci.no_gfactors)
                get_gfactors = false;
            else if(specification->ci.gfactors)
                get_gfactors = true;
            else
            {   if(nrID || num_solutions > 50)
                    get_gfactors = false;
                else
                    get_gfactors = true;
            }
        }

        if(get_gfactors)
        {
            GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);
            g_factors.CalculateGFactors(levelvec);
            levels->Store(hID, levelvec);
        }
    }

    // Set up output options
    bool ShowgFactors = true;
    if(specification->ci.no_gfactors)
    {   ShowgFactors = false;
    }

    bool ShowPercentages = true;
    if(specification->ci.output_no_configs || specification->ci.single_configuration_ci)
    {   ShowPercentages = false;
    }

    bool ShowRelConfigPercentages = specification->ci.output_print_relativistic_configurations;

    if(specification->ci.output_print_inline)
    {
        std::string sep = specification->ci.output_separator.value_or(" ");

        if(specification->ci.output_max_displayed_energy)
        {   // Truncate display at max energy
            double max_energy = specification->ci.output_max_displayed_energy.value();
            if(ShowRelConfigPercentages)
                levelvec.PrintInline<RelativisticConfiguration>(ShowPercentages, ShowgFactors, sep, max_energy);
            else
                levelvec.PrintInline(ShowPercentages, ShowgFactors, sep, max_energy);
        }
        else
        {   if(ShowRelConfigPercentages)
                levelvec.PrintInline<RelativisticConfiguration>(ShowPercentages, ShowgFactors, sep);
            else
                levelvec.PrintInline(ShowPercentages, ShowgFactors, sep);
        }
    }
    else
    {   double min_percent_displayed = 101.;
        if(ShowPercentages)
            min_percent_displayed = specification->ci.output_min_displayed_percent.value_or(1.0);

        if(specification->ci.output_max_displayed_energy)
        {   // Truncate display at max energy
            double DavidsonMaxEnergy = specification->ci.output_max_displayed_energy.value();
            if(ShowRelConfigPercentages)
                levelvec.Print<RelativisticConfiguration>(min_percent_displayed, DavidsonMaxEnergy);
            else
                levelvec.Print(min_percent_displayed, DavidsonMaxEnergy);
        }
        else
        {   if(ShowRelConfigPercentages)
                levelvec.Print<RelativisticConfiguration>(min_percent_displayed);
            else
                levelvec.Print(min_percent_displayed);
        }
    }

    return levelvec;
}

LevelVector Atom::SingleElectronConfigurations(pHamiltonianID sym)
{
    LevelVector levelvec = levels->GetLevels(sym);

    if(levelvec.NumLevels())
        return levelvec;

    if(hf_electron == nullptr)
        MakeCIIntegrals();

    OrbitalInfo info(dynamic_cast<SingleOrbitalID*>(sym.get())->GetOrbitalInfo());
    // Make relativistic configuration
    RelativisticConfiguration config;
    config.insert(std::make_pair(info, 1));
    config.GetProjections(angular_library, sym->GetSymmetry(), sym->GetTwoJ());

    angular_library->GenerateCSFs();

    // Make "level"
    levelvec.hID = sym;
    levelvec.configs = std::make_shared<RelativisticConfigList>(config);
    double energy = hf_electron->GetMatrixElement(info, info);
    levelvec.Resize(1, 1);
    levelvec.eigenvalues[0] = energy;
    levelvec.eigenvectors(0, 0) = 1;
    levels->Store(sym, levelvec);

    return levelvec;
}
}
