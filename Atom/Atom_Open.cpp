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

void Atom::RunOpen(bool include_mbpt)
{
    check_size_only = false;
    generate_mbpt_integrals = false;

    DebugOptions.LogFirstBuild(false);
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);
    DebugOptions.LogMBPT(false);
    //DebugOptions.LogAugerRate(true);

    core->ToggleOpenShellCore();

    //CreateCustomBasis(include_mbpt);
    //CreateRBasis(include_mbpt);
    CreateBSplineBasis(include_mbpt);

    DebugOptions.OutputHFExcited(false);

    if(include_mbpt)
        ReadOrWriteBasis();

    core->ToggleClosedShellCore();

    if(include_mbpt && generate_mbpt_integrals)
    {
        double mbpt_delta = 0.0;

        // Uncomment to include sigma3.
        //sigma3 = new Sigma3Calculator(lattice, core, excited);
        //sigma3->SetEnergyShift(mbpt_delta/Constant::HartreeEnergy_cm);

        //GenerateIntegralsMBPT(true, false, mbpt_delta);
        if(!check_size_only)
            CollateIntegralsMBPT(16);
    }
    else
        GenerateIntegrals(include_mbpt);

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
    core->SetNuclearInverseMass(0.);

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
    if(!integrals)
    {   integrals = new CIIntegralsMBPT(*excited);
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
        integrals->WriteOneElectronIntegrals();
        integrals->WriteTwoElectronIntegrals();
    }

    core->ToggleClosedShellCore();

    #ifdef _MPI
        // Wait for root node to finish writing, then update integrals.
        MPI::COMM_WORLD.Barrier();
    #endif

    integrals->Update();
}

void Atom::GenerateIntegrals(bool MBPT_CI)
{
    core->ToggleOpenShellCore();
    core->SetNuclearInverseMass(0.);
    if(MBPT_CI)
      Read();
    core->ToggleClosedShellCore();

    if(MBPT_CI)
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
    unsigned int two_j;

    // Even parity
    for(two_j =0; two_j <= 12; two_j += 2)
    {
        SymEigenstates[Symmetry(two_j, even)] = NULL;
    }

    // Odd parity
    for(two_j = 0; two_j <= 12; two_j += 2)
    {
        SymEigenstates[Symmetry(two_j, odd)] = NULL;
    }
}

ConfigGenerator* Atom::GenerateConfigurations(const Symmetry& sym)
{
    // Number of electron excitations (e.g. 2 for SD-CI) .
    unsigned int electron_excitations = 2;

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

        Configuration config1;
        config1.SetOccupancy(NonRelInfo(6, 2), 2);
        leading_configs.insert(config1);
/*
        Configuration config2;
        config2.SetOccupancy(NonRelInfo(3, 2), 2);
        config2.SetOccupancy(NonRelInfo(4, 1), 1);
        leading_configs.insert(config2);
/*
        Configuration config3;
        config3.SetOccupancy(NonRelInfo(3, 2), 2);
        config3.SetOccupancy(NonRelInfo(4, 1), 2);
        leading_configs.insert(config3);
/*
        Configuration config4;
        config4.SetOccupancy(NonRelInfo(2, 1), 1);
        config4.SetOccupancy(NonRelInfo(4, 3), 1);
        leading_configs.insert(config4);
*/
        generator->Clear();
        generator->AddLeadingConfigurations(leading_configs);
        
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
