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

namespace Ambit
{
void Atom::MakeMBPTIntegrals()
{
    bool make_new_integrals = !user_input.search(1, "--no-new-mbpt");
    bool check_sizes = user_input.search("--check-sizes");
    bool one_body_mbpt = user_input.search(3, "-s1", "-s12", "-s123");
    bool two_body_mbpt = user_input.search(4, "-s2", "-s12", "-s23", "-s123");

    // First make any Brueckner orbitals, and use these everywhere
    if(user_input.search("MBPT/--brueckner") && !check_sizes)
    {
        GenerateBruecknerOrbitals(make_new_integrals);
    }

    if(!make_new_integrals && !check_sizes)
        return;

    // Bare integrals for MBPT
    pHFIntegrals bare_one_body_integrals = std::make_shared<HFIntegrals>(orbitals, hf);
    pSlaterIntegrals bare_two_body_integrals = std::make_shared<SlaterIntegralsUnorderedDense>(orbitals, hartreeY);

    // MBPT calculators
    std::string fermi_orbitals = user_input("MBPT/EnergyDenomOrbitals", "");
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
    if(user_input.search("MBPT/--use-subtraction") || user_input.vector_variable_size("Basis/InjectOrbitals"))
        is_open_shell = true;
    if(user_input.search("MBPT/--no-subtraction"))
        is_open_shell = false;

    bool use_box = !user_input.search("MBPT/--no-extra-box");
    bool include_core = !user_input.search("MBPT/--no-core");
    mbpt_integrals_one->IncludeCore(include_core, include_core && is_open_shell);
    mbpt_integrals_two->IncludeCore(include_core, include_core && is_open_shell, include_core && use_box);

    bool include_valence = user_input.search("MBPT/--use-valence");
    mbpt_integrals_one->IncludeValence(include_valence && is_open_shell);
    mbpt_integrals_two->IncludeValence(include_valence, include_valence && is_open_shell, include_valence && use_box);

    // Adjust delta
    double delta = user_input("MBPT/Delta", 0.0);
    core_mbpt->SetEnergyShift(delta);
    val_mbpt->SetEnergyShift(delta);

    // Also set a floor for the valence-MBPT energy denominators
    double energy_denom_floor = user_input("MBPT/EnergyDenomFloor", 0.01);
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
        unsigned int num_limits = user_input.vector_variable_size("MBPT/TwoBody/StorageLimits");
        if(num_limits)
        {
            for(unsigned int i = 0; i < mmin(num_limits, 4); i++)
            {
                int max_pqn = user_input("MBPT/TwoBody/StorageLimits", 0, i);

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
            mbpt_integrals_one->Read(identifier + ".one.int");
            unsigned int size = mbpt_integrals_one->CalculateOneElectronIntegrals(valence, valence);
            *logstream << "One-body MBPT integrals complete: size = " << size << std::endl;
        }

        if(two_body_mbpt)
        {
            mbpt_integrals_two->Read(identifier + ".two.int");
            unsigned int size = mbpt_integrals_two->CalculateTwoElectronIntegrals(valence_subset[0], valence_subset[1], valence_subset[2], valence_subset[3]);
            *logstream << "Two-body MBPT integrals complete: size = " << size << std::endl;
        }
    }
}

void Atom::MakeIntegrals()
{
    ClearIntegrals();

    pSlaterIntegrals two_body_integrals;
    bool one_body_mbpt = user_input.search(3, "-s1", "-s12", "-s123");
    bool two_body_mbpt = user_input.search(4, "-s2", "-s12", "-s23", "-s123");
    bool three_body_mbpt = user_input.search(3, "-s3", "-s23", "-s123");

    if(!two_body_mbpt)
        two_body_integrals.reset(new SlaterIntegralsUnorderedDense(orbitals, hartreeY));
    else
        two_body_integrals.reset(new SlaterIntegralsUnorderedDense(orbitals, hartreeY, false));

    if(three_body_mbpt)
    {
        std::string fermi_orbitals = user_input("MBPT/EnergyDenomOrbitals", "");
        threebody_electron = std::make_shared<Sigma3Calculator>(orbitals, two_body_integrals, fermi_orbitals);
        threebody_electron->IncludeCore(!user_input.search("MBPT/--no-core"));
        threebody_electron->IncludeValence(user_input.search("MBPT/--use-valence"));
    }

    auto& valence = orbitals->valence;

    // We don't need the two body integrals if we are just checking sizes and they are not needed by ConfigGenerator.
    if(user_input.search("--check-sizes")
       && !user_input.VariableExists("CI/ConfigurationAverageEnergyRange")
       && !user_input.VariableExists("CI/SmallSide/ConfigurationAverageEnergyRange")
       && !user_input.search(2, "CI/--print-relativistic-configurations", "CI/--print-configurations")
       && !user_input.search(2, "CI/SmallSide/--print-relativistic-configurations", "CI/SmallSide/--print-configurations"))
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
        if(!user_input.search(2, "--no-ci", "--no-CI"))
        {
            two_body_integrals->CalculateTwoElectronIntegrals(valence, valence, valence, valence);
            if(user_input.search("--check-sizes"))
                *outstream << "\nNum Coulomb integrals: " << two_body_integrals->size() << std::endl;
        }

        // Add stored MBPT integrals
        if(one_body_mbpt)
        {
            unsigned int scaling_length = user_input.vector_variable_size("MBPT/OneBody/Scaling");
            if(scaling_length)
            {
                std::map<int, double> scaling;
                for(int i = 0; i < scaling_length-1; i+=2)
                {
                    int kappa = user_input("MBPT/OneBody/Scaling", 0, i);
                    double lambda = user_input("MBPT/OneBody/Scaling", 0.0, i+1);

                    scaling[kappa] = lambda;
                }

                hf_electron->Read(identifier + ".one.int", scaling);
            }
            else
                hf_electron->Read(identifier + ".one.int");
        }
        if(two_body_mbpt)
            two_body_integrals->Read(identifier + ".two.int");

        bool include_box = two_body_mbpt && !user_input.search("MBPT/--no-extra-box");
        bool include_off_parity = include_box || hartreeY->OffParityExists();
        twobody_electron = std::make_shared<TwoElectronCoulombOperator>(two_body_integrals, include_off_parity);

        // Calculate integrals for sigma3
        if(three_body_mbpt)
        {
            threebody_electron->UpdateIntegrals();
            if(user_input.search("--check-sizes"))
                *outstream << "\nSigma3 Coulomb integrals: " << threebody_electron->GetStorageSize() << std::endl;
        }
    }
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
        if(user_input.VariableExists("AngularDataDirectory"))
            angular_directory = user_input("AngularDataDirectory", "");

        angular_library = std::make_shared<AngularDataLibrary>(angular_directory);
    }
}

pLevelStore Atom::ChooseHamiltoniansAndRead(pAngularDataLibrary angular_lib)
{
    // Read existing levels?
    bool use_read = true;
    if(user_input.search(2, "--clean", "-c"))
        use_read = false;

    if(user_input.search("--configuration-average"))
    {
        ConfigGenerator gen(orbitals, user_input);
        if(twobody_electron == nullptr)
            MakeIntegrals();
        allconfigs = gen.GenerateConfigurations(hf_electron, twobody_electron->GetIntegrals());

        // Create empty level set
        levels = std::make_shared<LevelMap>(identifier, nullptr);
        return levels;
    }

    // Get angular library
    InitialiseAngularDataLibrary(angular_lib);

    if(user_input.search("CI/--memory-saver"))
    {
        std::string dir = user_input("LevelDirectory", "");
        if(dir.empty())
            levels = std::make_shared<FileSystemLevelStore>(identifier, angular_library);
        else
            levels = std::make_shared<FileSystemLevelStore>(dir, identifier, angular_library);
    }
    else if(use_read)
    {
        pHamiltonianID key;
        if(user_input.search(2, "--no-ci", "--no-CI"))
            key = std::make_shared<SingleOrbitalID>();
        else if(user_input.search(2, "CI/--single-configuration-ci", "CI/--single-configuration-CI"))
            key = std::make_shared<NonRelID>();
        else
            key = std::make_shared<HamiltonianID>();

        levels = std::make_shared<LevelMap>(key, identifier, angular_library);
    }
    else
        levels = std::make_shared<LevelMap>(identifier, angular_library);

    if(user_input.search(2, "--no-ci", "--no-CI"))
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
        ConfigGenerator gen(orbitals, user_input);

        if(twobody_electron == nullptr)
            MakeIntegrals();

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
    std::vector<int> even_symmetries;
    std::vector<int> odd_symmetries;
    pHamiltonianID key;

    // Get user symmetries if CI/--all-symmetries is not used
    if(!user_input.search("CI/--all-symmetries"))
    {
        // Even parity
        even_symmetries.resize(user_input.vector_variable_size("CI/EvenParityTwoJ"), 0);
        for(int i = 0; i < even_symmetries.size(); i++)
            even_symmetries[i] = user_input("CI/EvenParityTwoJ", 0, i);

        // Odd parity
        odd_symmetries.resize(user_input.vector_variable_size("CI/OddParityTwoJ"), 0);
        for(int i = 0; i < odd_symmetries.size(); i++)
            odd_symmetries[i] = user_input("CI/OddParityTwoJ", 0, i);
    }

    // Single configuration CI: each Hamiltonian only uses one non-rel configuration
    if(user_input.search(2, "CI/--single-configuration-ci", "CI/--single-configuration-CI"))
    {
        ConfigGenerator gen(orbitals, user_input);
        pConfigList nrlist = gen.GenerateNonRelConfigurations(rlist);

        if(user_input.search("CI/--all-symmetries"))
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
                int maxTwoJ = nrconfig.GetTwiceMaxProjection();

                if(P == Parity::even)
                {
                    for(int& two_j: even_symmetries)
                        if(two_j <= maxTwoJ)
                        {
                            key = std::make_shared<NonRelID>(nrconfig, two_j);
                            levels->keys.insert(key);
                        }
                }
                else
                {   for(int& two_j: odd_symmetries)
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
        if(user_input.search("CI/--all-symmetries"))
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

            for(int twoJ = mmax(max_twoJ_even%2, 0); twoJ <= max_twoJ_even; twoJ+=2)
                even_symmetries.push_back(twoJ);
            for(int twoJ = mmax(max_twoJ_odd%2, 0); twoJ <= max_twoJ_odd; twoJ+=2)
                odd_symmetries.push_back(twoJ);
        }

        // Populate levels
        for(int& two_j: even_symmetries)
        {   key = std::make_shared<HamiltonianID>(two_j, Parity::even);
            levels->keys.insert(key);
        }
        for(int& two_j: odd_symmetries)
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
    if(user_input.search(2, "--no-ci", "--no-CI"))
        return;

    // Get angular library
    InitialiseAngularDataLibrary(angular_lib);

    // Don't read existing levels
    levels = std::make_shared<LevelMap>(identifier, angular_library);

    // CI integrals
    MakeIntegrals();

    // Generate configurations again; don't read from disk. */
    ConfigGenerator gen(orbitals, user_input);
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
    if(user_input.search(2, "CI/--single-configuration-ci", "CI/--single-configuration-CI"))
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
    // This function is public and can call the other CalculateEnergies variants.
    LevelVector levelvec = levels->GetLevels(hID);
    std::string filename = identifier + ".levels";

    auto nrID = std::dynamic_pointer_cast<NonRelID>(hID);
    auto soID = std::dynamic_pointer_cast<SingleOrbitalID>(hID);

    if(soID)
    {
        if(levelvec.levels.empty())
            levelvec = SingleElectronConfigurations(hID);
    }
    else
    {   // Get relativistic configurations
        pRelativisticConfigList& configs = levelvec.configs;

        if(configs == nullptr)
        {
            ConfigGenerator gen(orbitals, user_input);

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
                        if(user_input.VariableExists("CI/ConfigurationAverageEnergyRange")
                           || user_input.VariableExists("CI/SmallSide/ConfigurationAverageEnergyRange")
                           || user_input.search(2, "CI/--print-relativistic-configurations", "CI/--print-configurations")
                           || user_input.search(2, "CI/SmallSide/--print-relativistic-configurations", "CI/SmallSide/--print-configurations"))
                        {
                            if(twobody_electron == nullptr)
                                MakeIntegrals();

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
        int num_solutions = user_input("CI/NumSolutions", 6);
        num_solutions = (num_solutions? mmin(num_solutions, configs->NumCSFs()): configs->NumCSFs());
        bool do_CI = (levelvec.levels.size() < num_solutions);
        if(do_CI)
        {
            if(twobody_electron == nullptr)
                MakeIntegrals();

            std::unique_ptr<HamiltonianMatrix> H;
            if(threebody_electron)
                H.reset(new HamiltonianMatrix(hf_electron, twobody_electron, threebody_electron, leading_configs, configs));
            else
                H.reset(new HamiltonianMatrix(hf_electron, twobody_electron, configs));

            // If we're using OpenMP then the chunksize should be a multiple of the number of threads
            int default_chunksize = 4;

            H->GenerateMatrix(user_input("CI/ChunkSize", default_chunksize));
            //H->PollMatrix();

            if(user_input.search("CI/Output/--write-hamiltonian"))
            {
                std::string hamiltonian_filename = identifier + "." + hID->Name() + ".matrix";

                // Convert spaces to underscores in filename
                std::replace_if(hamiltonian_filename.begin(), hamiltonian_filename.end(),
                                [](char c){ return (c =='\r' || c =='\t' || c == ' ' || c == '\n');}, '_');
                H->Write(hamiltonian_filename);
            }

            if(user_input.search("CI/Output/--print-hamiltonian"))
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

            #ifdef AMBIT_USE_SCALAPACK
            if(user_input.search("CI/--scalapack"))
            {
                if(user_input.VariableExists("CI/MaxEnergy"))
                {
                    double max_energy = user_input("CI/MaxEnergy", 0.0);
                    levelvec = H->SolveMatrixScalapack(hID, max_energy);
                }
                else
                {
                    levelvec = H->SolveMatrixScalapack(hID, num_solutions, false);
                }
            }
            else
            #endif
            levelvec = H->SolveMatrix(hID, num_solutions);
            levels->Store(hID, levelvec);
        }

        // Check if gfactor overrides are present, otherwise decide on course of action
        bool get_gfactors = levels->GFactorsNeeded();

        // Don't bother doing any checks if we've got pre-calculated g-factors stored
        if(get_gfactors)
        {
            if(hID->GetTwoJ() == 0 || user_input.search("CI/--no-gfactors"))
                get_gfactors = false;
            else if(user_input.search("CI/--gfactors"))
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
    if(user_input.search("CI/--no-gfactors"))
    {   ShowgFactors = false;
    }

    bool ShowPercentages = true;
    if(user_input.search("CI/Output/--no-configs") || user_input.search(2, "CI/--single-configuration-ci", "CI/--single-configuration-CI"))
    {   ShowPercentages = false;
    }

    bool ShowRelConfigPercentages = user_input.search("CI/Output/--print-relativistic-configurations");

    if(user_input.search("CI/Output/--print-inline"))
    {
        std::string sep = user_input("CI/Output/Separator", " ");

        if(user_input.search("CI/Output/MaxDisplayedEnergy"))
        {   // Truncate display at max energy
            double max_energy = user_input("CI/Output/MaxDisplayedEnergy", 0.);
            if(ShowRelConfigPercentages)
                levelvec.PrintInline<RelativisticConfiguration>(max_energy, ShowPercentages, ShowgFactors, sep);
            else
                levelvec.PrintInline(max_energy, ShowPercentages, ShowgFactors, sep);
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
            min_percent_displayed = user_input("CI/Output/MinimumDisplayedPercentage", 1.);

        if(user_input.search("CI/Output/MaxDisplayedEnergy"))
        {   // Truncate display at max energy
            double DavidsonMaxEnergy = user_input("CI/Output/MaxDisplayedEnergy", 0.);
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

    if(levelvec.levels.size())
        return levelvec;

    if(hf_electron == nullptr)
        MakeIntegrals();

    double eigenvector = 1.;

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
    levelvec.levels.emplace_back(std::make_shared<Level>(energy, &eigenvector, sym, 1));
    levels->Store(sym, levelvec);

    return levelvec;
}
}
