#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Atom.h"
#include "Universal/Constant.h"
#include "HartreeFock/NonRelInfo.h"
#include "Configuration/ConfigGenerator.h"
#include "Configuration/HamiltonianMatrix.h"
#include "Configuration/MPIHamiltonianMatrix.h"
#include "Basis/BSplineBasis.h"
#include "Basis/ReadBasis.h"

void Atom::RunMultipleElectron()
{
    check_size_only = userInput_.search("--check-size-only");
    generate_mbpt_integrals = userInput_.search("--generate-integrals-mbpt")
                             && userInput_.search("MBPTBasis");

    bool collate_mbpt_integrals = generate_mbpt_integrals || userInput_.search("--collate-integrals-mbpt");

    DebugOptions.LogFirstBuild(false);
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);
    DebugOptions.LogMBPT(false);
    //DebugOptions.LogAugerRate(true);

    core->ToggleClosedShellCore();

    if(!check_size_only && collate_mbpt_integrals)
        CollateIntegralsMBPT(NumProcessors);

    double mbpt_delta = userInput_("MBPTDelta", 0.0);

    // Uncomment to include sigma3.
    //sigma3 = new Sigma3Calculator(lattice, core, excited);
    //sigma3->SetEnergyShift(mbpt_delta/Constant::HartreeEnergy_cm);

    if(generate_mbpt_integrals)
    {   GenerateIntegralsMBPT(true, false, mbpt_delta);

        // Collate new integrals
        if(!check_size_only)
            CollateIntegralsMBPT(NumProcessors);
    }

    GenerateIntegrals();
    ChooseSymmetries();

    if(check_size_only)
        CheckMatrixSizes();
    else
    {   // Warning: Need to have generated integrals already.
        CalculateEnergies();
    }
}

void Atom::GenerateIntegralsMBPT(bool CoreMBPT, bool ValenceMBPT, double delta)
{
    core->ToggleOpenShellCore();

    integrals = new CIIntegralsMBPT(*excited);
    integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);

    Read();
    core->ToggleClosedShellCore();

    if(mbpt)
        delete mbpt;

    if(CoreMBPT)
    {   mbpt = new CoreMBPTCalculator(lattice, core, excited_mbpt);
        integralsMBPT->IncludeMBPT1(true, mbpt);
        integralsMBPT->IncludeMBPT2(true, mbpt);
        integralsMBPT->IncludeExtraBoxDiagrams(true);
    }
    else
        mbpt = NULL;

    if(valence_mbpt)
        delete valence_mbpt;
    
    if(ValenceMBPT)
    {   valence_mbpt = new ValenceCalculator(lattice, core, excited_mbpt);
        integralsMBPT->IncludeValenceMBPT1(true, valence_mbpt);
        integralsMBPT->IncludeValenceMBPT2(true, valence_mbpt);
        integralsMBPT->IncludeValenceExtraBoxDiagrams(true);
    }
    else
        valence_mbpt = NULL;

    integralsMBPT->SetTwoElectronStorageLimits(7, 7);

    // Affects both core and valence MBPT if extra box diagrams are included.
    // To include box diagrams in Hamiltonian, uncomment the #defines at the top of HamiltonianMatrix.cpp.
    integralsMBPT->SetExtraBoxDiagramLimits(7, 7);

    if(mbpt)
        mbpt->SetEnergyShift(delta/Constant::HartreeEnergy_cm);
    if(valence_mbpt)
        valence_mbpt->SetEnergyShift(delta/Constant::HartreeEnergy_cm);

    integrals->IncludeValenceSMS(false);
    integrals->SetIdentifier(identifier);    
    integrals->Clear();
    if(!check_size_only)
        integrals->Update();
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
    {   integralsMBPT->ReadMultipleOneElectronIntegrals(identifier, num_processors);
        integralsMBPT->ReadMultipleTwoElectronIntegrals(identifier, num_processors);
        integrals->WriteOneElectronIntegrals(true);
        integrals->WriteTwoElectronIntegrals(true);
    }

    core->ToggleClosedShellCore();

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

    core->ToggleOpenShellCore();
    core->SetNuclearInverseMass(0.);
    //if(MBPT_CI)
    //  Read();
    core->ToggleClosedShellCore();

    // Are we using MBPT?
    if(userInput_.search("--use-mbpt"))
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
    num_symmetries = userInput_.vector_variable_size("EvenParityTwoJ");
    for(i = 0; i < num_symmetries; i++)
    {
        two_j = userInput_("EvenParityTwoJ", 0, i);
        SymEigenstates[Symmetry(two_j, even)] = NULL;
    }

    // Odd parity
    num_symmetries = userInput_.vector_variable_size("OddParityTwoJ");
    for(i = 0; i < num_symmetries; i++)
    {
        two_j = userInput_("OddParityTwoJ", 0, i);
        SymEigenstates[Symmetry(two_j, odd)] = NULL;
    }

    if(SymEigenstates.empty())
    {   *errstream << "USAGE: No symmetries requested (EvenParityTwoJ or OddParityTwoJ)" << std::endl;
        exit(1);
    }
}

ConfigGenerator* Atom::GenerateConfigurations(const Symmetry& sym)
{
    // Number of electron excitations (e.g. 2 for SD-CI) .
    unsigned int electron_excitations = userInput_("ElectronExcitations", 2);

    // Generate non-relativistic configs from file.
    bool GenerateFromFile = false;
    
    ConfigGenerator* generator;

    if(GenerateFromFile)
        generator = new ConfigFileGenerator(excited, identifier, sym);
    else
        generator = new ConfigGenerator(excited, identifier, sym);

    bool read_from_disk = false;
    if(save_configurations)
        read_from_disk = generator->Read();

    if(!read_from_disk)
    {
        std::set<Configuration> leading_configs;

        int num_configs = userInput_.vector_variable_size("LeadingConfigurations");

        if(num_configs < 1)
        {   *errstream << "USAGE: Need LeadingConfigurations (e.g. 3d7, 4s2 3d5)" << std::endl;
            exit(1);
        }

        generator->Clear();
        for(int i = 0; i < num_configs; i++)
        {
            const std::string name = userInput_("LeadingConfigurations", "", i);
            Configuration config(name);

            // Check that the configuration gels with the number of electrons
            if(config.NumParticles() != numValenceElectrons_)
            {   *errstream << "USAGE: LeadingConfiguration " << name
                           << " does not have correct number of valence electrons." << std::endl;
                exit(1);
            }
            generator->AddLeadingConfiguration(config);
        }

        if(GenerateFromFile)
        {   ConfigFileGenerator* filegenerator = dynamic_cast<ConfigFileGenerator*>(generator);
            filegenerator->SetInputFile("PercentagesIn.txt");
            filegenerator->ReadConfigs(0.05);
        }
        else
        {   generator->GenerateMultipleExcitationsFromLeadingConfigs(electron_excitations);
        }

        generator->GenerateRelativisticConfigs();
        generator->GenerateProjections();

        if(save_configurations)
            generator->Write();
    }
    
    return generator;
}

void Atom::CheckMatrixSizes()
{
    // Two electron integral storage size
    if(integrals)
        *outstream << "\nNum coulomb integrals: " << integrals->GetStorageSize() << std::endl;

    SymmetryEigenstatesMap::iterator it = SymEigenstates.begin();

    while(it != SymEigenstates.end())
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
    SymmetryEigenstatesMap::iterator it = SymEigenstates.begin();

    while(it != SymEigenstates.end())
    {
        ConfigGenerator* conf_gen = GenerateConfigurations(it->first);
        Eigenstates* E = new Eigenstates(identifier, conf_gen);

        if(!E->Read())
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

            #ifdef _SCALAPACK
                H->WriteToFile("temp.matrix");
                MPIHamiltonianMatrix* MpiH = dynamic_cast<MPIHamiltonianMatrix*>(H);
                MpiH->SolveScalapack("temp.matrix", -1.47, *E, true);
            #else
                H->SolveMatrix(NumSolutions, *E, true);
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

        it->second = E;
        if(save_eigenstates)
            E->Clear();
        it++;
    }
}

Eigenstates* Atom::GetEigenstates(const Symmetry& sym)
{
    SymmetryEigenstatesMap::iterator it = SymEigenstates.find(sym);    
    if(it != SymEigenstates.end())
        return it->second;
    else
        return NULL;
}
