#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Atom.h"
#include "Universal/Constant.h"
#include "Universal/Enums.h"
#include "HartreeFock/NonRelInfo.h"
#include "Configuration/ConfigGenerator.h"
#include "Configuration/HamiltonianMatrix.h"
#include "Configuration/MPIHamiltonianMatrix.h"
#include "Configuration/MPIMatrix.h"
#include "Basis/BSplineBasis.h"
#include "Basis/ReadBasis.h"

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

    //core->ToggleOpenShellCore();
    //core->SetNuclearInverseMass(0.);
    //if(MBPT_CI)
    //  Read();
    core->ToggleClosedShellCore();

    if(includeSigma1 || includeSigma2)
        integrals = new CIIntegralsMBPT(*excited);
    else
        integrals = new CIIntegrals(*excited);

    integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);

    integrals->IncludeValenceSMS(false);
    integrals->SetIdentifier(identifier);
    integrals->Clear();
    if(!check_size_only)
    {   integrals->Update();

        if(sigma3)
            sigma3->UpdateIntegrals(excited);
    }
}

void Atom::ChooseSymmetries()
{
    int i, num_symmetries;
    int two_j;

    // Even parity
    num_symmetries = userInput_.vector_variable_size("CI/EvenParityTwoJ");
    for(i = 0; i < num_symmetries; i++)
    {
        two_j = userInput_("CI/EvenParityTwoJ", 0, i);
        symEigenstates[Symmetry(two_j, even)] = NULL;
    }

    // Odd parity
    num_symmetries = userInput_.vector_variable_size("CI/OddParityTwoJ");
    for(i = 0; i < num_symmetries; i++)
    {
        two_j = userInput_("CI/OddParityTwoJ", 0, i);
        symEigenstates[Symmetry(two_j, odd)] = NULL;
    }

    if(symEigenstates.empty())
    {   *errstream << "USAGE: No symmetries requested (EvenParityTwoJ or OddParityTwoJ)" << std::endl;
        exit(1);
    }
}

ConfigGenerator* Atom::GenerateConfigurations(const Symmetry& sym)
{
    bool allow_different_excitations = false;
    unsigned int electron_excitations = 0;
    int num_electron_excitation_inputs = userInput_.vector_variable_size("CI/ElectronExcitations");

    // Default if no input detected is 2 electron excitations.
    // If input is detected it can either be a number, which would set electron_excitations and use ValenceBasis to determine states to excite to,
    // otherwise if input is of the form 'X, Y, ...' where X is the pqn and Y is the string with the basis (eg. 8spdf) then the code will
    // allow X excitations to the states in Y
    if(num_electron_excitation_inputs < 1) 
    {
        electron_excitations = 2;
    }
    else if(num_electron_excitation_inputs == 1 && ((int) userInput_("CI/ElectronExcitations", 2)) >= 0)
    {
        electron_excitations = (int) userInput_("CI/ElectronExcitations", 2);
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

    // Generate non-relativistic configs from file.
    bool GenerateFromFile = false;

    ConfigGenerator* generator;

    if(GenerateFromFile)
        generator = new ConfigFileGenerator(excited, original_id, sym);
    else
        generator = new ConfigGenerator(excited, original_id, sym);

    bool read_from_disk = false;
    if(useRead)
        read_from_disk = generator->Read();

    if(!read_from_disk)
    {
        std::set<Configuration> leading_configs;

        int num_configs = userInput_.vector_variable_size("CI/LeadingConfigurations");

        if(num_configs < 1)
        {   *errstream << "USAGE: Need CI/LeadingConfigurations (e.g. '3d7, 4s2 3d5')" << std::endl;
            exit(1);
        }

        generator->Clear();
        for(int i = 0; i < num_configs; i++)
        {
            const std::string name = userInput_("CI/LeadingConfigurations", "", i);
            Configuration config(name);

            // Check that the configuration gels with the number of electrons
            if(config.NumParticles() != numValenceElectrons_)
            {   *errstream << "USAGE: LeadingConfiguration " << name
                           << " does not have correct number of valence electrons." << std::endl;
                exit(1);
            }
            generator->AddLeadingConfiguration(config);
        }

        // Adds extra configurations not to be considered leading configurations.
        if(userInput_.vector_variable_size("CI/ExtraConfigurations") > 0)
        {
            int num_extra_configs = userInput_.vector_variable_size("CI/ExtraConfigurations");
            for(int i = 0; i < num_extra_configs; i++) 
            {
                const std::string extraname = userInput_("CI/ExtraConfigurations", "", i);
                Configuration extraconfig(extraname);

                if(extraconfig.NumParticles() != numValenceElectrons_) 
                {
                    *errstream << "USAGE: LeadingConfiguration " << extraname
                           << " does not have correct number of valence electrons." << std::endl;
                    exit(1);
                }
                generator->AddNonRelConfiguration(extraconfig);
            }
        }

        if(GenerateFromFile)
        {   ConfigFileGenerator* filegenerator = dynamic_cast<ConfigFileGenerator*>(generator);
            filegenerator->SetInputFile("PercentagesIn.txt");
            filegenerator->ReadConfigs(0.05);
        }
        else if(allow_different_excitations)
        {
            unsigned int CI_num_excitations;
            std::vector<unsigned int> CI_electron_excitation_states;
            std::string CI_basis_string;

            NonRelInfoSet nrset;
            ConstStateIterator si = core->GetConstStateIterator();

            for(int i = 0; i < num_electron_excitation_inputs; i += 2)
            {
                CI_num_excitations = atoi(userInput_("CI/ElectronExcitations", "", i).c_str());
                CI_basis_string = userInput_("CI/ElectronExcitations", "", i+1);

                if(!ParseBasisSize(CI_basis_string.c_str(), CI_electron_excitation_states))
                {
                    *errstream << "USAGE: CI/ElectronExcitations = " << CI_basis_string << " incorrectly specified." << std::endl;
                    exit(1);
                }

                nrset.clear();
                nrset.AddConfigs(CI_basis_string.c_str());
                // Erase core set
                for(si.First(); !si.AtEnd(); si.Next())
                {
                    nrset.erase(si.GetOrbitalInfo());
                }

                generator->GenerateMultipleExcitationsFromLeadingConfigs(CI_num_excitations, &nrset);
            }
        }
        else
        {
            generator->GenerateMultipleExcitationsFromLeadingConfigs(electron_excitations);
        }

        generator->GenerateRelativisticConfigs();
        generator->GenerateProjections();

        if(useWrite)
            generator->Write();
    }

    return generator;
}

void Atom::CheckMatrixSizes()
{
    // Two electron integral storage size
    if(integrals)
        *outstream << "\nNum coulomb integrals: " << integrals->GetStorageSize() << std::endl;

    SymmetryEigenstatesMap::iterator it = symEigenstates.begin();

    while(it != symEigenstates.end())
    {
        // Generate configurations again; don't read from disk. */
        ConfigGenerator* conf_gen = GenerateConfigurations(it->first);
        Eigenstates* E = new Eigenstates(identifier, conf_gen);

        *outstream << "\nJ = " << it->first.GetJ() << ", P = ";
        if(it->first.GetParity() == even)
            *outstream << "even" << std::endl;
        else
            *outstream << "odd" << std::endl;

        // Print number of relativistic configurations
        *outstream << " Number of rel configurations = "
                   << E->GetConfigGenerator()->GetRelConfigs()->size() << std::endl;
        
        // Print number of JStates
        *outstream << " Number of J-configurations = "
                   << E->GetConfigGenerator()->GetNumJStates() << std::endl;

        delete E;
        it++;
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

Eigenstates* Atom::GetEigenstates(const Symmetry& sym)
{
    SymmetryEigenstatesMap::iterator it = symEigenstates.find(sym);    
    if(it != symEigenstates.end())
        return it->second;
    else
        return NULL;
}
