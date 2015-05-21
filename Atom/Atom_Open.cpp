#ifdef _MPI
#include <mpi.h>
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

void Atom::MakeMBPTIntegrals()
{
    bool make_new_integrals = user_input.search(2, "-m", "--mbpt");
    bool check_sizes = user_input.search("--check-sizes");
    bool one_body_mbpt = user_input.search(3, "-s1", "-s12", "-s123");
    bool two_body_mbpt = user_input.search(3, "-s2", "-s12", "-s123");

    // First make any Brueckner orbitals, and use these everywhere
    if(user_input.search("--brueckner") && !check_sizes)
    {
        pBruecknerDecorator brueckner(new BruecknerDecorator(hf));
        bool use_fg = user_input.search("MBPT/--brueckner-use-lower");
        bool use_gg = user_input.search("MBPT/--brueckner-use-lower-lower");
        brueckner->IncludeLower(use_fg, use_gg);

        // Attempt to read all requested kappas
        std::set<int> valence_kappas;
        for(auto& pair: *orbitals->valence)
            valence_kappas.insert(pair.first.Kappa());

        for(int kappa: valence_kappas)
            brueckner->Read(identifier, kappa);

        // Replace hf operator for rest of calculation
        hf = brueckner;

        // Make new sigma potentials if they haven't been read (slowly)
        if(make_new_integrals)
        {
            for(int kappa: valence_kappas)
            {   brueckner->CalculateSigma(kappa, orbitals, hartreeY);
                brueckner->Write(identifier, kappa);
            }
        }

        // And finally change all our valence orbitals to Brueckner orbitals
        pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
        pODESolver ode_solver(new AdamsSolver(integrator));
        HartreeFocker hartree_focker(ode_solver);

        DebugOptions.LogHFIterations(true);
        DebugOptions.OutputHFExcited(true);

        pOrbital brueckner_orbital;
        for(auto& pair: *orbitals->valence)
        {
            // Make a copy of the old orbital
            brueckner_orbital.reset(new Orbital(*pair.second));

            // Iterate
            hartree_focker.SolveOrbital(brueckner_orbital, brueckner);

            // Copy back to orbital manager
            *pair.second = *brueckner_orbital;
        }
    }

    if(!make_new_integrals && !check_sizes)
        return;

    // Bare integrals for MBPT
    pOneElectronIntegrals bare_one_body_integrals(new OneElectronIntegrals(orbitals, hf));
    pSlaterIntegrals bare_two_body_integrals(new SlaterIntegralsMap(orbitals, hartreeY));
    pCoreMBPTCalculator core_mbpt(new CoreMBPTCalculator(orbitals, bare_one_body_integrals, bare_two_body_integrals));

    auto& valence = orbitals->valence;

    pOneElectronMBPT mbpt_integrals_one(new OneElectronMBPT(orbitals, core_mbpt, hf));
    pCoreValenceIntegralsMap mbpt_integrals_two(new CoreValenceIntegralsMap(orbitals, core_mbpt));

    // Use subtraction diagrams and extra box diagrams?
    // First find out whether the core is our closed shell core, while we're at it, get maximum pqn.
    int max_pqn_in_core = 0;
    bool is_open_shell = false;
    for(auto pair: *orbitals->core)
    {
        double occ = open_core->GetOccupancy(pair.first);
        if(fabs(occ - double(pair.first.MaxNumElectrons())) > 0.01)
            is_open_shell = true;

        max_pqn_in_core = mmax(max_pqn_in_core, pair.first.PQN());
    }

    bool use_box = !user_input.search("MBPT/--no-box-diagrams");
    mbpt_integrals_one->IncludeCore(true, is_open_shell);
    mbpt_integrals_two->IncludeCore(true, is_open_shell, use_box);

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
        unsigned int num_limits = user_input.vector_variable_size("MBPT/TwoElectronStorageLimits");
        if(num_limits)
        {
            for(unsigned int i = 0; i < mmin(num_limits, 4); i++)
            {
                int max_pqn = user_input("MBPT/TwoElectronStorageLimits", 0, i);

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
        *outstream << "Num stored coulomb integrals for MBPT = " << core_mbpt->GetStorageSize();

        unsigned int size = mbpt_integrals_one->CalculateOneElectronIntegrals(valence, valence, true);
        *outstream << "\nNum one-body mbpt integrals: " << size;

        size = mbpt_integrals_two->CalculateTwoElectronIntegrals(valence_subset[0], valence_subset[1], valence_subset[2], valence_subset[3], true);
        *outstream << "\nNum two-body mbpt integrals: " << size << std::endl;
    }
    else
    {
        if(one_body_mbpt)
        {
            mbpt_integrals_one->CalculateOneElectronIntegrals(valence, valence);
            mbpt_integrals_one->Write(identifier + ".one.int");
        }

        if(two_body_mbpt)
        {
            mbpt_integrals_two->CalculateTwoElectronIntegrals(valence_subset[0], valence_subset[1], valence_subset[2], valence_subset[3]);
            mbpt_integrals_two->Write(identifier + ".two.int");
        }
    }
}

void Atom::MakeIntegrals()
{
    ClearIntegrals();

    pSlaterIntegrals two_body_integrals;
    bool one_body_mbpt = user_input.search(3, "-s1", "-s12", "-s123");
    bool two_body_mbpt = user_input.search(3, "-s2", "-s12", "-s123");

    if(!two_body_mbpt)
        two_body_integrals.reset(new SlaterIntegralsMap(orbitals, hartreeY));
    else
        two_body_integrals.reset(new SlaterIntegralsMap(orbitals, hartreeY, false));

    auto& valence = orbitals->valence;

    if(user_input.search("--check-sizes"))
    {
        unsigned int size = two_body_integrals->CalculateTwoElectronIntegrals(valence, valence, valence, valence, true);
        *outstream << "\nNum coulomb integrals: " << size << std::endl;
    }
    else
    {   // Normal integrals
        hf_electron.reset(new OneElectronIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(valence, valence);

        // Don't need two body integrals if we're not doing CI
        if(!user_input.search(2, "--no-ci", "--no-CI"))
            two_body_integrals->CalculateTwoElectronIntegrals(valence, valence, valence, valence);

        // Add stored MBPT integrals
        if(one_body_mbpt)
            hf_electron->Read(identifier + ".one.int");
        if(two_body_mbpt)
            two_body_integrals->Read(identifier + ".two.int");

        twobody_electron.reset(new TwoElectronCoulombOperator<pSlaterIntegrals>(two_body_integrals));
    }
}

void Atom::ClearIntegrals()
{
    hf_electron = nullptr;
    twobody_electron = nullptr;
}

void Atom::InitialiseAngularDataLibrary(pAngularDataLibrary trial)
{
    int num_valence_electrons;
    if(user_input.search("NumValenceElectrons"))
        num_valence_electrons = user_input("NumValenceElectrons", 0);
    else if(user_input.search(2, "--no-ci", "--no-CI"))
        num_valence_electrons = 1;
    else if(user_input.vector_variable_size("CI/LeadingConfigurations") > 0)
    {
        NonRelConfiguration nrconfig(user_input("CI/LeadingConfigurations", "", 0));
        num_valence_electrons = nrconfig.ElectronNumber();
    }

    if(trial && trial->GetElectronNumber() == num_valence_electrons)
        angular_library = trial;
    else if(angular_library == nullptr || angular_library->GetElectronNumber() != num_valence_electrons)
    {
        std::string angular_directory = string_macro(ANGULAR_DATA_DIRECTORY);
        if(user_input.search("AngularDataDirectory"))
            angular_directory = user_input("AngularDataDirectory", "");

        angular_library = std::make_shared<AngularDataLibrary>(num_valence_electrons, angular_directory);
    }
}

pLevelMap Atom::ChooseHamiltoniansAndRead(pAngularDataLibrary angular_lib)
{
    // Read existing levels?
    bool use_read = true;
    if(user_input.search(2, "--clean", "-c"))
        use_read = false;

    // Get angular library
    InitialiseAngularDataLibrary(angular_lib);

    if(use_read)
    {
        std::string filename = identifier + ".levels";

        pHamiltonianID key;
        if(user_input.search(2, "--no-ci", "--no-CI"))
            key = std::make_shared<SingleOrbitalID>();
        else if(user_input.search(2, "CI/--single-configuration-ci", "CI/--single-configuration-CI"))
            key = std::make_shared<NonRelID>();
        else
            key = std::make_shared<HamiltonianID>();

        levels = ReadLevelMap(key, filename, angular_library);
    }
    else
        levels = std::make_shared<LevelMap>();

    if(user_input.search(2, "--no-ci", "--no-CI"))
    {
        // Use all symmetries from valence set
        pHamiltonianID key;
        for(auto& it: *orbitals->valence)
        {
            key = std::make_shared<SingleOrbitalID>(it.first);
            (*levels)[key];
        }
    }
    else
    {   // Generate non-rel configurations and hence choose Hamiltonians
        ConfigGenerator gen(orbitals, user_input);
        nrconfigs = gen.GenerateNonRelConfigurations(hf, hartreeY);

        ChooseHamiltonians(nrconfigs);
    }

    return levels;
}

pLevelMap Atom::ChooseHamiltonians(pConfigList nrlist)
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
        if(user_input.search("CI/--all-symmetries"))
        {
            // Populate levels with all symmetries
            for(auto& nrconfig: *nrlist)
            {
                int maxTwoJ = nrconfig.GetTwiceMaxProjection();
                for(int two_j = maxTwoJ%2; two_j <= maxTwoJ; two_j += 2)
                {
                    key = std::make_shared<NonRelID>(nrconfig, two_j);
                    (*levels)[key];
                }
            }
        }
        else
        {   // Populate levels with symmetries found in sets.
            for(auto& nrconfig: *nrlist)
            {
                Parity P = nrconfig.GetParity();
                int maxTwoJ = nrconfig.GetTwiceMaxProjection();
                
                if(P == Parity::even)
                {
                    for(int& two_j: even_symmetries)
                        if(two_j <= maxTwoJ)
                        {
                            key = std::make_shared<NonRelID>(nrconfig, two_j);
                            (*levels)[key];
                        }
                }
                else
                {   for(int& two_j: odd_symmetries)
                        if(two_j <= maxTwoJ)
                        {
                            key = std::make_shared<NonRelID>(nrconfig, two_j);
                            (*levels)[key];
                        }
                }

            }
            
        }
    }
    // Standard CI: all non-rel configs together
    else
    {   // Get symmetries
        if(user_input.search("CI/--all-symmetries"))
        {
            int max_twoJ_even = -1;
            int max_twoJ_odd = -1;

            for(auto& nrconfig: *nrlist)
            {
                if(nrconfig.GetParity() == Parity::even)
                    max_twoJ_even = mmax(max_twoJ_even, nrconfig.GetTwiceMaxProjection());
                else
                    max_twoJ_odd = mmax(max_twoJ_odd, nrconfig.GetTwiceMaxProjection());
            }

            for(int twoJ = mmax(max_twoJ_even%2, 0); twoJ <= max_twoJ_even; twoJ+=2)
                even_symmetries.push_back(twoJ);
            for(int twoJ = mmax(max_twoJ_odd%2, 0); twoJ <= max_twoJ_odd; twoJ+=2)
                odd_symmetries.push_back(twoJ);
        }

        // Populate levels
        for(int& two_j: even_symmetries)
        {   key = std::make_shared<HamiltonianID>(two_j, Parity::even);
            (*levels)[key];
        }
        for(int& two_j: odd_symmetries)
        {   key = std::make_shared<HamiltonianID>(two_j, Parity::odd);
            (*levels)[key];
        }
    }

    if(levels->empty())
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

    levels = std::make_shared<LevelMap>();

    // Generate configurations again; don't read from disk. */
    ConfigGenerator gen(orbitals, user_input);
    nrconfigs = gen.GenerateNonRelConfigurations(hf, hartreeY);

    // CI integrals
    MakeIntegrals();
    ChooseHamiltonians(nrconfigs);

    // Get complete list of relativistic configs for all symmetries
    std::map<Symmetry, pRelativisticConfigList> all_relconfigs;
    for(auto& pair: *levels)
    {
        pHamiltonianID key = pair.first;
        Symmetry sym = pair.first->GetSymmetry();
        pRelativisticConfigList configs;

        if(NonRelID* nrid = dynamic_cast<NonRelID*>(key.get()))
        {
            pConfigList nrconfiglist = std::make_shared<ConfigList>(nrid->GetNonRelConfiguration());
            configs = gen.GenerateRelativisticConfigurations(nrconfiglist, nrid->GetSymmetry());
        }
        else
        {
            configs = gen.GenerateRelativisticConfigurations(nrconfigs, key->GetSymmetry());
        }

        key->SetRelativisticConfigList(configs);

        pRelativisticConfigList& symrel = all_relconfigs[sym];
        if(!symrel)
            symrel = std::make_shared<RelativisticConfigList>(*configs);
        else
            symrel->append(*configs);
    }

    // Generate CSFs for each symmetry
    *outstream << "Calculating CSFs: " << std::endl;
    unsigned int total_levels = 0;
    for(auto& pair: all_relconfigs)
    {
        const Symmetry& sym = pair.first;
        *outstream << "J(P) = " << sym.GetJ() << "(" << ShortName(sym.GetParity()) << "): "
                   << std::setw(6) << std::right << pair.second->size() << " rel. configurations; "
                   << std::flush;

        // Generate all projections for this symmetry and write
        gen.GenerateProjections(pair.second, sym, sym.GetTwoJ(), angular_library);

        *outstream << pair.second->NumCSFs() <<  " CSFs." << std::endl;
        total_levels += pair.second->NumCSFs();
    }
    *outstream << "\nTotal number of levels (all symmetries included) = " << total_levels << std::endl;

    // Get Hamiltonian sizes
    if(user_input.search(2, "CI/--single-configuration-ci", "CI/--single-configuration-CI"))
    {
        *outstream << "\nHamiltonian matrix sizes: " << std::endl;
        for(auto& pair: *levels)
        {
            auto rconfigs = pair.first->GetRelativisticConfigList();
            unsigned int num_CSFs = 0;
            for(auto& config: *rconfigs)
            {
                num_CSFs += angular_library->GetData(config, pair.first->GetSymmetry(), pair.first->GetTwoJ())->NumCSFs();
            }

            *outstream << pair.first->Print() << ":"
                       << std::setw(6) << std::right << rconfigs->size() << " rel. configurations; "
                       << num_CSFs <<  " CSFs." << std::endl;
        }
    }
}

pLevelMap Atom::CalculateEnergies()
{
    ChooseHamiltoniansAndRead();

    for(auto& pair: *levels)
        CalculateEnergies(pair.first);

    return levels;
}

const LevelVector& Atom::CalculateEnergies(pHamiltonianID hID)
{
    // This function is public and can call the other CalculateEnergies variants.
    LevelVector& levelvec = (*levels)[hID];
    pHamiltonianID key = levels->find(hID)->first;

    std::string filename = identifier + ".levels";

    if(dynamic_cast<SingleOrbitalID*>(key.get()))
    {
        if(levelvec.size() == 0)
            SingleElectronConfigurations(key);
    }
    else
    {   // Get relativistic configurations
        pRelativisticConfigList configs = key->GetRelativisticConfigList();

        if(configs == nullptr)
        {   configs = hID->GetRelativisticConfigList();
            key->SetRelativisticConfigList(configs);
        }

        if(configs == nullptr)
        {
            ConfigGenerator gen(orbitals, user_input);

            if(NonRelID* nrid = dynamic_cast<NonRelID*>(key.get()))
            {
                pConfigList nrconfiglist = std::make_shared<ConfigList>(nrid->GetNonRelConfiguration());
                configs = gen.GenerateRelativisticConfigurations(nrconfiglist, nrid->GetSymmetry(), angular_library);
            }
            else
            {
                if(nrconfigs == nullptr)
                    nrconfigs = gen.GenerateNonRelConfigurations(hf, hartreeY);
                configs = gen.GenerateRelativisticConfigurations(nrconfigs, key->GetSymmetry(), angular_library);
            }

            key->SetRelativisticConfigList(configs);
        }

        // Only continue if we don't have enough levels
        int num_solutions = NumSolutions? mmin(NumSolutions, configs->NumCSFs()): configs->NumCSFs();
        if(levelvec.size() < num_solutions)
        {
            if(twobody_electron == nullptr)
                MakeIntegrals();

            HamiltonianMatrix H(hf_electron, twobody_electron, configs);

//        if(sigma3)
//            H->IncludeSigma3(sigma3);

            H.GenerateMatrix();
            //H->PollMatrix();

            if(user_input.search("--write-hamiltonian"))
            {
                std::string hamiltonian_filename = identifier + "." + key->Name() + ".matrix";

                // Convert spaces to underscores in filename
                std::replace_if(hamiltonian_filename.begin(), hamiltonian_filename.end(),
                                [](char c){ return (c =='\r' || c =='\t' || c == ' ' || c == '\n');}, '_');
                H.Write(hamiltonian_filename);
            }

            if((user_input("CI/Output/PrintH", "false") == "true") || (user_input("CI/Output/PrintH", 0) == 1))
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
                *outstream << "Matrix Before:\n" << H << std::endl;
            }

            #ifdef _SCALAPACK
                H->WriteToFile("temp.matrix");
                MPIHamiltonianMatrix* MpiH = dynamic_cast<MPIHamiltonianMatrix*>(H);
                MpiH->SolveScalapack("temp.matrix", MaxEnergy, *E, true, num_solutions);
            #else
                levelvec = H.SolveMatrix(key, num_solutions);
            #endif

            if(key->GetTwoJ() != 0)
            {
                GFactorCalculator g_factors(hf->GetOPIntegrator(), orbitals);
                g_factors.CalculateGFactors(levelvec);
            }

//            ConfigFileGenerator* filegenerator = dynamic_cast<ConfigFileGenerator*>(conf_gen);
//            if(filegenerator)
//            {   filegenerator->SetOutputFile("PercentagesOut.txt");
//                filegenerator->WriteConfigs();
//            }

            WriteLevelMap(*levels, filename);
        }
    }

    // Set up output options
    bool ShowgFactors = true;
    if((user_input("CI/Output/ShowgFactors", "true") == "false") || (user_input("CI/Output/ShowgFactors", 1) == 0))
    {   ShowgFactors = false;
    }

    bool ShowPercentages = true;
    if((user_input("CI/Output/ShowPercentages", "true") == "false") || (user_input("CI/Output/ShowPercentages", 1) == 0))
    {   ShowPercentages = false;
    }
    double min_percent_displayed;

    if(ShowPercentages)
    {   min_percent_displayed = user_input("CI/Output/MinimumDisplayedPercentage", 1.);
    }
    else
    {   min_percent_displayed = 101.;
    }

    if(user_input.search("CI/Output/MaxDisplayedEnergy"))
    {   // Truncate display at max energy
        double DavidsonMaxEnergy = user_input("CI/Output/MaxDisplayedEnergy", 0.);
        Print(levelvec, min_percent_displayed, DavidsonMaxEnergy);
    }
    else
        Print(levelvec, min_percent_displayed);

    return levelvec;
}

const LevelVector& Atom::SingleElectronConfigurations(pHamiltonianID sym)
{
    LevelVector& levelvec = (*levels)[sym];

    if(levelvec.size())
        return levelvec;

    if(hf_electron == nullptr)
        MakeIntegrals();

    double eigenvector = 1.;

    OrbitalInfo info(dynamic_cast<SingleOrbitalID*>(sym.get())->GetOrbitalInfo());
    // Make relativistic configuration
    RelativisticConfiguration config;
    config.insert(std::make_pair(info, 1));
    config.GetProjections(angular_library, sym->GetSymmetry(), sym->GetTwoJ());

    // Make "level"
    pRelativisticConfigList configlist(new RelativisticConfigList(config));
    sym->SetRelativisticConfigList(configlist);
    double energy = hf_electron->GetMatrixElement(info, info);
    pLevel level(std::make_shared<Level>(energy, &eigenvector, sym, 1));

    levelvec = LevelVector(1, level);

    return levelvec;
}

/*
void Atom::RunMultipleElectron()
{
    check_size_only = userInput_.search("--check-sizes");
    generate_mbpt_integrals = userInput_.search("--generate-integrals-mbpt")
                             && (mbptBasisString != "");

    // Don't generate MBPT integrals if we're not saving them!
    if(generate_mbpt_integrals && !useWrite && !check_size_only)
    {   *outstream << "USAGE: Cannot use \"--generate-integrals-mbpt\" with \"--dont-save\" unless\n"
                   << "       only doing a \"--check-sizes\"." << std::endl;
        exit(1);
    }

    bool collate_mbpt_integrals = (generate_mbpt_integrals || userInput_.search("--collate-integrals-mbpt"))
                                  && !check_size_only;

    if(collate_mbpt_integrals)
    {
        while(!RunIndexAtEnd())
        {   CollateIntegralsMBPT(NumProcessors);
            RunIndexNext(false);
        }
        RunIndexBegin(false);
    }

    // MBPT Options
    bool includeSall = userInput_.search(2, "-s123", "--include-sigma-all");
    includeSigma1 = userInput_.search(3, "-s1", "-s12", "--include-sigma1");
    includeSigma2 = userInput_.search(3, "-s2", "-s12", "--include-sigma2");
    includeSigma3 = userInput_.search(2, "-s3", "--include-sigma3");

    includeSigma1 = includeSigma1 || includeSall;
    includeSigma2 = includeSigma2 || includeSall;
    includeSigma3 = includeSigma3 || includeSall;

    // If not specified for mbpt integral calculation, make both Sigma 1 and Sigma 2
    if(generate_mbpt_integrals && !includeSigma1 && !includeSigma2)
    {   includeSigma1 = true;
        includeSigma2 = true;
    }

    unsigned int num_mbpt_delta = userInput_.vector_variable_size("MBPT/Delta");
    if((num_mbpt_delta > 1) && (num_mbpt_delta != multiple_length))
    {   *outstream << "ERROR: MBPT/Delta has too few elements. Must match multiple run length." << std::endl;
        exit(1);
    }
    for(unsigned int i = 0; i < num_mbpt_delta; i++)
    {   mbpt_delta.push_back(userInput_("MBPT/Delta", 0.0, i));
    }

    // Generate the MBPT integrals
    if(generate_mbpt_integrals)
    {
        InitialiseIntegralsMBPT(true, false);

        if(check_size_only)
        {   *outstream << std::endl;
            unsigned int stored_size = integrals->GetStorageSize();
            *outstream << "Total calculated MBPT integrals: " << stored_size << std::endl;
        }
        else
        {   // Generate all MBPT integrals
            while(!RunIndexAtEnd())
            {
                SetRunParameters(false);
                SetRunCore();
                SetRunIntegrals();  // This does the actual calculations

                RunIndexNext(false);
            }
            RunIndexBegin(false);

            // Collate new integrals
            if(collate_mbpt_integrals)
            {
                while(!RunIndexAtEnd())
                {   CollateIntegralsMBPT(NumProcessors);
                    RunIndexNext(false);
                }
                RunIndexBegin(false);
            }
        }
    }

    // We have finished generating MBPT integrals and collated them.
    // Create the integral sets for creating Hamiltonian
    if(integrals)
        delete integrals;

    if(includeSigma2)
        integrals = new CIIntegralsMBPT(*excited);
    else
        integrals = new CIIntegrals(*excited);
    integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);

    if(includeSigma3)
        sigma3 = new Sigma3Calculator(lattice, core, excited);

    SetRunParameters(false);
    SetRunCore();
    SetRunIntegrals(true);

    // Choose J, Parity symmetries for CI run.
    ChooseSymmetries();

    if(check_size_only)
        CheckMatrixSizes();
    else
    {   // Warning: Need to have generated integrals already.
        #ifdef _SCALAPACK
            MaxEnergy = userInput_("CI/MaxEnergy", 0.0);
            if(userInput_.search("CI/MaxEnergy"))
                NumSolutions = userInput_("CI/NumSolutions", 0);
            else
                NumSolutions = userInput_("CI/NumSolutions", 6);
        #else
            NumSolutions = userInput_("CI/NumSolutions", 6);
        #endif

        CalculateEnergies();
    }
}

void Atom::InitialiseIntegralsMBPT(bool CoreMBPT, bool ValenceMBPT)
{
    integrals = new CIIntegralsMBPT(*excited);
    integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);

    core->ToggleClosedShellCore();

    if(mbpt)
        delete mbpt;

    if(CoreMBPT)
    {   mbpt = new CoreMBPTCalculator(lattice, core, excited_mbpt);
        integralsMBPT->IncludeMBPT1(includeSigma1, mbpt);
        integralsMBPT->IncludeMBPT2(includeSigma2, mbpt);
        integralsMBPT->IncludeExtraBoxDiagrams(includeSigma2);
    }
    else
        mbpt = NULL;

    if(valence_mbpt)
        delete valence_mbpt;
    
    if(ValenceMBPT)
    {   valence_mbpt = new ValenceCalculator(lattice, core, excited_mbpt);
        integralsMBPT->IncludeValenceMBPT1(includeSigma1, valence_mbpt);
        integralsMBPT->IncludeValenceMBPT2(includeSigma2, valence_mbpt);
        integralsMBPT->IncludeValenceExtraBoxDiagrams(includeSigma2);
    }
    else
        valence_mbpt = NULL;

    std::vector<int> two_electron_limits;
    for(unsigned int i = 0; i < 3; i++)
        two_electron_limits.push_back(userInput_("MBPT/TwoElectronStorageLimits", 0, i));
    integralsMBPT->SetTwoElectronStorageLimits(two_electron_limits[0], two_electron_limits[1], two_electron_limits[2]);

    // Affects both core and valence MBPT if extra box diagrams are included.
    // To include box diagrams in Hamiltonian, uncomment the #defines at the top of HamiltonianMatrix.cpp.
    if(userInput_.vector_variable_size("MBPT/BoxDiagramStorageLimits") >= 1)
    {   two_electron_limits.clear();
        for(unsigned int i = 0; i < 3; i++)
            two_electron_limits.push_back(userInput_("MBPT/BoxDiagramStorageLimits", 0, i));
    }
    integralsMBPT->SetExtraBoxDiagramLimits(two_electron_limits[0], two_electron_limits[1], two_electron_limits[2]);
}

void Atom::CollateIntegralsMBPT(unsigned int num_processors)
{
    bool integrals_previously_exist = true;

    if(!integrals)
    {   integrals_previously_exist = false;
        integrals = new CIIntegralsMBPT(*excited);
        integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);
    }
    else if(integralsMBPT)
    {   // Stop doing MBPT calculations
        integralsMBPT->IncludeMBPT1(false);
        integralsMBPT->IncludeMBPT2(false);
        integralsMBPT->IncludeExtraBoxDiagrams(false);
        integralsMBPT->IncludeValenceMBPT1(false);
        integralsMBPT->IncludeValenceMBPT2(false);
        integralsMBPT->IncludeValenceExtraBoxDiagrams(false);
    }

    integrals->SetIdentifier(identifier);
    integrals->Clear();

    if(ProcessorRank == 0)
    {   if(includeSigma1)
        {   integralsMBPT->ReadMultipleOneElectronIntegrals(identifier, num_processors);
            integrals->WriteOneElectronIntegrals(true);
        }
        if(includeSigma2)
        {   integralsMBPT->ReadMultipleTwoElectronIntegrals(identifier, num_processors);
            integrals->WriteTwoElectronIntegrals(true);
        }
    }

    #ifdef _MPI
        // Wait for root node to finish writing
        MPI::COMM_WORLD.Barrier();
    #endif

    if(!integrals_previously_exist)
    {   delete integrals;
        integrals = NULL;
        integralsMBPT = NULL;
    }
}

void Atom::GenerateIntegrals()
{
    if(integrals)
        delete integrals;

    if(includeSigma1 || includeSigma2)
        integrals = new CIIntegralsMBPT(*excited);
    else
        integrals = new CIIntegrals(*excited);

    integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);

    integrals->IncludeValenceSMS(false);

    integrals = new CIIntegrals(hf, hartreeY, orbitals->valence, identifier, hartreeY_reverse_symmetry);

    if(!user_input.search("--check-sizes"))
    {
        integrals->Update();
        if(sigma3)
            sigma3->UpdateIntegrals(excited);
    }
}

void Atom::CalculateEnergies()
{
    // Set up output options
    bool ShowgFactors = true;
    if((userInput_("CI/Output/ShowgFactors", "true") == "false") || (userInput_("CI/Output/ShowgFactors", 1) == 0))
    {   ShowgFactors = false;
    }

    bool ShowPercentages = true;
    if((userInput_("CI/Output/ShowPercentages", "true") == "false") || (userInput_("CI/Output/ShowPercentages", 1) == 0))
    {   ShowPercentages = false;
    }
    double min_percent_displayed;
    
    if(ShowPercentages)
    {   min_percent_displayed = userInput_("CI/Output/MinimumDisplayedPercentage", 1.);
    }
    else
    {   min_percent_displayed = 101.;
    }

    bool TruncateDisplayAtMaxEnergy = userInput_.search("CI/Output/MaxDisplayedEnergy");
    double DavidsonMaxEnergy = 0.;
    if(TruncateDisplayAtMaxEnergy)
    {
        DavidsonMaxEnergy = userInput_("CI/Output/MaxDisplayedEnergy", 0.);
    }

    // Create and solve Hamiltonian matrix for all symmetries
    SymmetryEigenstatesMap::iterator it = symEigenstates.begin();

    while(it != symEigenstates.end())
    {
        ConfigGenerator* conf_gen = GenerateConfigurations(it->first);

        bool give_conf_gen_to_eigenstates = false;
        if(NumberRunsSelected() == 1)
            give_conf_gen_to_eigenstates = true;

        while(!RunIndexAtEnd())
        {
            Eigenstates* E = new Eigenstates(identifier, conf_gen, give_conf_gen_to_eigenstates);

            SetRunParameters(true);
            SetRunCore();
            SetRunIntegrals();

            if(!useRead || !E->Read())
            {
                HamiltonianMatrix* H;

                #ifdef _MPI
                    H = new MPIHamiltonianMatrix(*integrals, conf_gen);
                #else
                    H = new HamiltonianMatrix(*integrals, conf_gen);
                #endif

                if(sigma3)
                    H->IncludeSigma3(sigma3);

                H->GenerateMatrix();
                //H->PollMatrix();

                if((userInput_("CI/Output/PrintH", "false") == "true") || (userInput_("CI/Output/PrintH", 0) == 1))
                {
                    #ifdef _MPI
                        std::string filename = identifier + "." + it->first.GetString() + ".matrix";
                        dynamic_cast<MPIHamiltonianMatrix*>(H)->WriteToFile(filename, false);
                    #else
                        RelativisticConfigList::iterator rel_it = conf_gen->GetRelConfigs()->begin();
                        while(rel_it != conf_gen->GetRelConfigs()->end())
                        {
                            *outstream << rel_it->Name();
                            if(rel_it++ != conf_gen->GetRelConfigs()->end())
                            {
                                *outstream << ",";
                            }
                        }
                        *outstream << std::endl;
                        
                        *outstream << std::setprecision(12);
                        *outstream << "Matrix Before:" << std::endl;
                        for(unsigned int i = 0; i < H->GetMatrix()->GetSize(); i++)
                        {
                            for(unsigned int j = 0; j < H->GetMatrix()->GetSize(); j++)
                            {
                                *outstream << H->GetMatrix()->At(i,j) << " ";
                            }
                            *outstream << std::endl;
                        }
                    #endif
                }
                
                #ifdef _SCALAPACK
                    H->WriteToFile("temp.matrix");
                    MPIHamiltonianMatrix* MpiH = dynamic_cast<MPIHamiltonianMatrix*>(H);
                    MpiH->SolveScalapack("temp.matrix", MaxEnergy, *E, true, NumSolutions);
                #else
                    H->SolveMatrix(NumSolutions, *E, GetSolutionMap(), ShowgFactors, TruncateDisplayAtMaxEnergy, min_percent_displayed, DavidsonMaxEnergy);
                #endif

                delete H;

                ConfigFileGenerator* filegenerator = dynamic_cast<ConfigFileGenerator*>(conf_gen);
                if(filegenerator)
                {   filegenerator->SetOutputFile("PercentagesOut.txt");
                    filegenerator->WriteConfigs();
                }

                if(save_eigenstates)
                    E->Write();
            }
            else
            {   E->Print();
            }

            // Keep E if not multiple run and save_eigenstates is true
            if(save_eigenstates && NumberRunsSelected() == 1)
            {   it->second = E;
                //E->Clear();
            }
            else
            {   //delete E;
            }

            RunIndexNext(false);
        }

        if(!give_conf_gen_to_eigenstates)
            delete conf_gen;

        RunIndexBegin(false);
        it++;
    }
}
*/
