#include "Include.h"
#include "Atom.h"
#include "Universal/Constant.h"
#include "HartreeFock/NonRelInfo.h"
#include "Configuration/ConfigGenerator.h"
#include "Configuration/HamiltonianMatrix.h"
#include "Configuration/MPIHamiltonianMatrix.h"
#include "Basis/BSplineBasis.h"

#ifdef _MPI
#include <mpi.h>
#endif

void Atom::RunOpen()
{
    DebugOptions.LogFirstBuild(false);
    DebugOptions.LogHFIterations(false);
    //DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);
    DebugOptions.LogMBPT(false);

    //CreateCustomBasis();
    //CreateRBasis();
    CreateBSplineBasis();

    DebugOptions.OutputHFExcited(false);

    double mbpt_delta = 0.0;

    GenerateIntegrals(false);
    ChooseSymmetries();

    // Uncomment to include sigma3.
    //sigma3 = new Sigma3Calculator(lattice, core, excited);
    //sigma3->SetEnergyShift(mbpt_delta/Constant::HartreeEnergy_cm);

    //CheckMatrixSizes();

    // Warning: Need to have generated integrals already.
    CalculateEnergies();
}

void Atom::GenerateIntegralsMBPT(bool CoreMBPT, bool ValenceMBPT, double delta)
{
    integrals = new CIIntegralsMBPT(*excited);
    integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);

    // Create excited state basis. Should be a superset of the CI basis.
    excited_mbpt = new BSplineBasis(lattice, core);
    excited_mbpt->SetIdentifier(&identifier);
    dynamic_cast<BSplineBasis*>(excited_mbpt)->SetParameters(40, 7, 25.);
    std::vector<unsigned int> num_states_per_l;
    num_states_per_l.push_back(40);
    num_states_per_l.push_back(40);
    num_states_per_l.push_back(39);
    num_states_per_l.push_back(38);
    num_states_per_l.push_back(37);
    excited_mbpt->CreateExcitedStates(num_states_per_l);

    if(CoreMBPT)
    {   mbpt = new MBPTCalculator(lattice, core, excited_mbpt);
        integralsMBPT->IncludeMBPT1(true, mbpt);
        integralsMBPT->IncludeMBPT2(true, mbpt);
        integralsMBPT->IncludeExtraBoxDiagrams(true);
    }
    else
    {   if(mbpt)
            delete mbpt;
        mbpt = NULL;
    }

    if(ValenceMBPT)
    {   valence_mbpt = new ValenceCalculator(lattice, core, excited_mbpt);
        integralsMBPT->IncludeValenceMBPT1(true, valence_mbpt);
        integralsMBPT->IncludeValenceMBPT2(true, valence_mbpt);
        integralsMBPT->IncludeValenceExtraBoxDiagrams(true);
    }
    else
    {   if(valence_mbpt)
            delete valence_mbpt;
        valence_mbpt = NULL;
    }

    // Affects both core and valence MBPT if extra box diagrams are included.
    // To include box diagrams in Hamiltonian, uncomment the #defines at the top of HamiltonianMatrix.cpp.
    integralsMBPT->SetExtraBoxDiagramLimits(4, 4);

    if(mbpt)
        mbpt->SetEnergyShift(delta/Constant::HartreeEnergy_cm);
    if(valence_mbpt)
        valence_mbpt->SetEnergyShift(delta/Constant::HartreeEnergy_cm);

    integrals->IncludeValenceSMS(false);
    integrals->SetIdentifier(identifier);    
    integrals->Clear();
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

    #ifdef _MPI
        // Wait for root node to finish writing, then update integrals.
        //MPI::Intracomm& comm_world = MPI::COMM_WORLD;
        //comm_world.Barrier();
        MPI::COMM_WORLD.Barrier();
        integrals->Update();
    #endif
}

void Atom::GenerateIntegrals(bool MBPT_CI)
{
    if(MBPT_CI)
        integrals = new CIIntegralsMBPT(*excited);
    else
        integrals = new CIIntegrals(*excited);

    integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);

    integrals->IncludeValenceSMS(false);
    integrals->SetIdentifier(identifier);
    integrals->Clear();
    integrals->Update();
}

void Atom::ChooseSymmetries()
{
    unsigned int two_j;

    // Even parity
    for(two_j = 0; two_j <= 2; two_j += 2)
    {
        SymEigenstates[Symmetry(two_j, even)] = NULL;
    }

    // Odd parity
    for(two_j = 0; two_j <= 2; two_j+=2)
    {
        SymEigenstates[Symmetry(two_j, odd)] = NULL;
    }
}

Eigenstates* Atom::GenerateConfigurations(const Symmetry& sym, bool try_read)
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
    if(try_read)
        read_from_disk = generator->Read();

    if(!read_from_disk)
    {
        std::set<Configuration> leading_configs;

        Configuration config1;
        config1.SetOccupancy(NonRelInfo(2, 0), 2);
        //config1.SetOccupancy(NonRelInfo(2, 1), 1);
        leading_configs.insert(config1);

        //Configuration config2;
        //config2.SetOccupancy(NonRelInfo(3, 0), 1);
        //config2.SetOccupancy(NonRelInfo(3, 1), 1);
        //leading_configs.insert(config2);

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
        
        generator->Write();
    }

    Eigenstates* E = new Eigenstates(identifier, generator);
    
    return E;
}

void Atom::CheckMatrixSizes()
{
    // Two electron integral storage size
    *outstream << "\nNum coulomb integrals: " << integrals->GetStorageSize() << std::endl;

    SymmetryEigenstatesMap::iterator it = SymEigenstates.begin();

    while(it != SymEigenstates.end())
    {
        // Generate configurations again; don't read from disk. */
        Eigenstates* E = GenerateConfigurations(it->first, false);

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
    core->ToggleOpenShellCore();
    core->SetNuclearInverseMass(0.);
    //Read();
    Write();
    core->ToggleClosedShellCore();

    SymmetryEigenstatesMap::iterator it = SymEigenstates.begin();

    while(it != SymEigenstates.end())
    {
        Eigenstates* E = GenerateConfigurations(it->first);

        if(!E->Read())
        {
            HamiltonianMatrix* H;

            #ifdef _MPI
                H = new MPIHamiltonianMatrix(*integrals, E->GetConfigGenerator());
            #else
                H = new HamiltonianMatrix(*integrals, E->GetConfigGenerator());
            #endif

            if(sigma3)
                H->IncludeSigma3(sigma3);

            H->GenerateMatrix();
            //H->PollMatrix();

            #ifdef _SCALAPACK
                H->WriteToFile("temp.matrix");
                MPIHamiltonianMatrix* MpiH = dynamic_cast<MPIHamiltonianMatrix*>(H);
                MpiH->SolveScalapack("temp.matrix", 0.5, *E, true);
            #else
                H->SolveMatrix(NumSolutions, *E, true);
            #endif

            delete H;

            ConfigFileGenerator* filegenerator = dynamic_cast<ConfigFileGenerator*>(E->GetConfigGenerator());
            if(filegenerator)
            {   filegenerator->SetOutputFile("PercentagesOut.txt");
                filegenerator->WriteConfigs();
            }

            E->Write();
        }

        it->second = E;        
        //E->Clear();
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

/*
void Atom::RunOpen()
{
    DebugOptions.LogFirstBuild(false);
    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);
    DebugOptions.LogMBPT(false);

    //CreateCustomBasis();
    //CreateRBasis();
    CreateBSplineBasis();
    DebugOptions.OutputHFExcited(false);

    SD_CI = true;
    MBPT_CI = false;

    if(MBPT_CI)
    {   integrals = new CIIntegralsMBPT(*excited);
        integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);

        excited_mbpt = new BSplineBasis(lattice, core);
        excited_mbpt->SetIdentifier(&identifier);
        dynamic_cast<BSplineBasis*>(excited_mbpt)->SetParameters(40, 7, 25.);
        std::vector<unsigned int> num_states_per_l;
        num_states_per_l.push_back(40);
        num_states_per_l.push_back(40);
        num_states_per_l.push_back(39);
        num_states_per_l.push_back(38);
        num_states_per_l.push_back(37);
        excited_mbpt->CreateExcitedStates(num_states_per_l);

        valence_mbpt = new ValenceCalculator(lattice, core, excited_mbpt);

        integralsMBPT->IncludeValenceMBPT1(true, valence_mbpt);
        integralsMBPT->IncludeValenceMBPT2(true, valence_mbpt);
        integralsMBPT->IncludeValenceExtraBoxDiagrams(true);//, 4, 4);
    }
    else
    {   integrals = new CIIntegrals(*excited);
        integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);
    }

    //integralsMBPT->SetTwoElectronStorageLimits(4, 4);
    //sigma3 = new Sigma3Calculator(lattice, core, excited);

    unsigned int NumParticles = 2;
//    CheckMatrixSizes();
//    return;

    unsigned int two_j;
    NumSolutions = 6;

    *outstream << "\nGS Parity:\n" << std::endl;
    for(two_j = 0; two_j <= 0; two_j += 2)
    {
        HamiltonianMatrix* H = CreateHamiltonian(two_j, NumParticles, even);
        //H->IncludeSigma3(sigma3);

        OpenShellEnergy(two_j, H);
        //DoOpenShellSMS(two_j, H);
        //DoOpenShellVolumeShift(two_j, H);

        delete H;
    }
    return;

    NumSolutions = 6;

    *outstream << "\nOpposite Parity:\n" << std::endl;
    for(two_j = 0; two_j <= 8; two_j+=2)
    {
        HamiltonianMatrix* H = CreateHamiltonian(two_j, NumParticles, odd);
        //H->IncludeSigma3(sigma3);
        
        OpenShellEnergy(two_j, H);
        //DoOpenShellSMS(two_j, H);
        //DoOpenShellVolumeShift(two_j, H);

        delete H;
    }
}

HamiltonianMatrix* Atom::CreateHamiltonian(unsigned int twoJ, unsigned int num_particles, Parity P)
{
    unsigned int electron_excitations;
    if(SD_CI)
        electron_excitations = 2;
    else
        electron_excitations = num_particles;

    std::set<Configuration> leading_configs;

    Configuration config1;
    config1.SetOccupancy(NonRelInfo(2, 0), 2);
    leading_configs.insert(config1);

    //Configuration config2;
    //config2.SetOccupancy(NonRelInfo(3, 0), 1);
    //config2.SetOccupancy(NonRelInfo(3, 1), 1);
    //leading_configs.insert(config2);

    GenerateFromFile = false;
    generator = new ConfigGenerator(excited, identifier, twoJ, P);

    generator->ClearConfigLists();
    generator->AddLeadingConfigurations(leading_configs);
    
    if(GenerateFromFile)
    {   ConfigFileGenerator* filegenerator = dynamic_cast<ConfigFileGenerator*>(generator);
        filegenerator->SetInputFile("PercentagesIn.txt");
        filegenerator->ReadConfigs(P, 0.05);
    }
    else
    {   generator->GenerateMultipleExcitationsFromLeadingConfigs(electron_excitations, P);
    }

    generator->GenerateRelativisticConfigs();
    generator->GenerateProjections(twoJ);

    generator->Write();

    HamiltonianMatrix* H;

    #ifdef _MPI
        H = new MPIHamiltonianMatrix(*integrals, generator);
    #else
        H = new HamiltonianMatrix(*integrals, generator);
    #endif

    return H;
}

void Atom::CheckMatrixSizes()
{
    // Two electron integral storage size
    *outstream << "Num coulomb integrals: " << integrals->GetStorageSize() << std::endl;

    unsigned int NumParticles = 2;
    unsigned int two_j;

    *outstream << "\nGS Parity:\n" << std::endl;
    for(two_j = 0; two_j <= 10; two_j += 2)
    {
        *outstream << "J = " << two_j/2. << std::endl;
        HamiltonianMatrix* H = CreateHamiltonian(two_j, NumParticles, even);
        *outstream << " Number of rel configurations = " << generator->GetRelConfigs()->size() << std::endl;
        delete H;
    }

    *outstream << "\nOpposite Parity:\n" << std::endl;
    for(two_j = 0; two_j <= 8; two_j+=2)
    {
        *outstream << "J = " << two_j/2. << std::endl;
        HamiltonianMatrix* H = CreateHamiltonian(two_j, NumParticles, odd);
        *outstream << " Number of rel configurations = " << generator->GetRelConfigs()->size() << std::endl;
        delete H;
    }
}

void Atom::OpenShellEnergy(int twoJ, HamiltonianMatrix* H)
{
    core->ToggleOpenShellCore();
    core->SetNuclearInverseMass(0.);
    //Read();
    //Write();
    core->ToggleClosedShellCore();

    integrals->IncludeValenceSMS(false);
    integrals->SetIdentifier(identifier);
    integrals->Clear();
    //integralsMBPT->ReadMultipleOneElectronIntegrals(identifier, 1);
    //integralsMBPT->ReadMultipleTwoElectronIntegrals(identifier, 1);
    //integrals->WriteOneElectronIntegrals();
    //integrals->WriteTwoElectronIntegrals();
    integrals->Update();

    H->GenerateMatrix();
    //H->WriteToFile("temp.matrix");
    //H->PollMatrix();

    #ifdef _SCALAPACK
        MPIHamiltonianMatrix* MpiH = dynamic_cast<MPIHamiltonianMatrix*>(H);
        MpiH->SolveScalapack("temp.matrix", 0.5, twoJ, true);
    #else
        H->SolveMatrix(NumSolutions, twoJ, true);
    #endif

    ConfigFileGenerator* filegenerator = dynamic_cast<ConfigFileGenerator*>(generator);
    if(filegenerator)
    {   filegenerator->SetOutputFile("PercentagesOut.txt");
        filegenerator->WriteConfigs();
    }
}

void Atom::DoOpenShellSMS(int twoJ, HamiltonianMatrix* H)
{
    integrals->IncludeValenceSMS(true);
    std::string original_id = identifier;

    // delta = E_CI - E_HF
    const static double delta[] = {-153032.342993, -153171.872643,
                                   -153311.988576, -153452.692228,
                                   -153593.985096};
    int p_delta;
    
    for(double ais = -0.002; ais <= 0.002; ais += 0.001)
    {
        *outstream << "\nNuclearInverseMass = " << ais << std::endl;

        std::stringstream ss;
        ss << (int)(ais * 1000);
        identifier = original_id + '_' + ss.str();
        
        core->ToggleOpenShellCore();
        core->SetNuclearInverseMass(ais);
        //core->Update();
        //excited->Update();
        //if(MBPT_CI)
        //    excited_mbpt->Update();
        //Write();
        Read();
        core->ToggleClosedShellCore();

        p_delta = (int)(ais * 1000) + 2;
        if(mbpt)
            mbpt->SetEnergyShift(delta[p_delta]/Constant::HartreeEnergy_cm);
        if(sigma3)
            sigma3->SetEnergyShift(delta[p_delta]/Constant::HartreeEnergy_cm);

        integrals->SetIdentifier(identifier);
        //integrals->Clear();
        //integralsMBPT->ReadMultipleTwoElectronIntegrals(identifier, 16);
        //integrals->WriteTwoElectronIntegrals();
        integrals->Update();

        H->GenerateMatrix();

        if(ais == 0.0)  // Get g-factors
            H->SolveMatrix(NumSolutions, twoJ, true);
        else
            H->SolveMatrix(NumSolutions, twoJ);
    }

    identifier = original_id;
}

void Atom::SMS_V0(int twoJ, HamiltonianMatrix* H)
{
    integrals->IncludeValenceSMS(false);
    for(double ais = -0.002; ais <= 0.002; ais += 0.004)
    {
        *outstream << "\nNuclearInverseMass = " << ais << std::endl;

        core->ToggleOpenShellCore();
        core->SetNuclearInverseMass(ais);
        core->Update();
        excited->Update();
        core->SetNuclearInverseMass(0.);
        core->ToggleClosedShellCore();

        integrals->Update();
        H->GenerateMatrix();
        H->SolveMatrix(NumSolutions, twoJ);
    }
}

void Atom::SMS_V1(int twoJ, HamiltonianMatrix* H)
{
    core->ToggleOpenShellCore();
    core->SetNuclearInverseMass(0.);
    core->Update();
    excited->Update();
    core->ToggleClosedShellCore();

    integrals->IncludeValenceSMS(false);
    for(double ais = -0.002; ais <= 0.002; ais += 0.004)
    {
        *outstream << "\nNuclearInverseMass = " << ais << std::endl;

        core->SetNuclearInverseMass(ais);

        integrals->Update();
        H->GenerateMatrix();
        H->SolveMatrix(NumSolutions, twoJ);
    }
}

void Atom::SMS_V2(int twoJ, HamiltonianMatrix* H)
{
    OpenShellEnergy(twoJ, H);

    // Just do one point since it is just a matrix element anyway.
    for(double ais = 0.002; ais <= 0.002; ais += 0.004)
    {
        *outstream << "\nNuclearInverseMass = " << ais << std::endl;

        core->SetNuclearInverseMass(ais);
        H->GetEigenvalues();
    }
}

void Atom::DoOpenShellFSModifyR(int twoJ, HamiltonianMatrix* H)
{
    integrals->IncludeValenceSMS(false);
    std::string original_id = identifier;

    // delta = E_CI - E_HF
    const static double delta[] = {-153309.654578, -153310.828119,
				   -153311.988571, -153313.136414,
				   -153314.272120};
    int p_delta;

    *outstream << "\nThickness = " << std::setprecision(3) << core->GetNuclearThickness()*Constant::AtomicToFermi;

    for(double r = 0.; r <= 12.; r += 3.)
    {
        *outstream << "\nRadius = " << std::setprecision(3) << r << std::endl;

        std::stringstream ss;
        ss << (int)(r);
        identifier = original_id + '_' + ss.str();

        core->ToggleOpenShellCore();
        core->SetNuclearRadius(r/Constant::AtomicToFermi);
        core->Update();
        excited->Update();
        //if(MBPT_CI)
        //    excited_mbpt->Update();
        //Write();
        //Read();
        core->ToggleClosedShellCore();

        p_delta = (int)(r/3. - 0.99);
        if(mbpt)
            mbpt->SetEnergyShift(delta[p_delta]/Constant::HartreeEnergy_cm);
        if(sigma3)
            sigma3->SetEnergyShift(delta[p_delta]/Constant::HartreeEnergy_cm);

        integrals->SetIdentifier(identifier);
        //integrals->Clear();
        //integralsMBPT->ReadMultipleOneElectronIntegrals(identifier, 16);
        //integrals->WriteOneElectronIntegrals();
        integrals->Update();

        H->GenerateMatrix();

        if(r == 3.0)  // Get g-factors
            H->SolveMatrix(NumSolutions, twoJ, true);
        else
            H->SolveMatrix(NumSolutions, twoJ);
    }

    identifier = original_id;
}

void Atom::DoOpenShellVolumeShift(int twoJ, HamiltonianMatrix* H)
{
    double deltaR = 0.05;
    static bool calculated_potential = false;
    if(!calculated_potential)
    {   core->CalculateVolumeShiftPotential(deltaR/Constant::AtomicToFermi);
        calculated_potential = true;
    }

    integrals->IncludeValenceSMS(false);
    std::string original_id = identifier;

    // delta = E_CI - E_HF
    const static double delta[] = {-153309.654578, -153310.828119,
				   -153311.988571, -153313.136414,
				   -153314.272120};
    int p_delta;

    *outstream << "\nThickness = " << std::setprecision(3) << core->GetNuclearThickness()*Constant::AtomicToFermi;
    *outstream << "\nRadius    = " << std::setprecision(3) << core->GetNuclearRadius()*Constant::AtomicToFermi;
    *outstream << "\nd(Radius) = " << std::setprecision(3) << deltaR;

    for(double ais = -100.; ais <= 100.; ais += 50.)
    {
        *outstream << "\nVolumeShiftParameter = " << std::setprecision(3) << ais << std::endl;

        std::stringstream ss;
        ss << (int)(ais);
        identifier = original_id + '_' + ss.str();

        core->ToggleOpenShellCore();
        core->SetVolumeShiftParameter(ais);
        //core->Update();
        //excited->Update();
        //if(MBPT_CI)
        //    excited_mbpt->Update();
        //Write();
        Read();
        core->ToggleClosedShellCore();

        p_delta = (int)(ais/50. + 2.01);
        if(mbpt)
            mbpt->SetEnergyShift(delta[p_delta]/Constant::HartreeEnergy_cm);
        if(sigma3)
            sigma3->SetEnergyShift(delta[p_delta]/Constant::HartreeEnergy_cm);

        integrals->SetIdentifier(identifier);
        //integrals->Clear();
        //integralsMBPT->ReadMultipleOneElectronIntegrals(identifier, 16);
        //integrals->WriteOneElectronIntegrals();
        integrals->Update();

        H->GenerateMatrix();

        if(ais == 0.0)  // Get g-factors
            H->SolveMatrix(NumSolutions, twoJ, true);
        else
            H->SolveMatrix(NumSolutions, twoJ);
    }

    identifier = original_id;
}

void Atom::DoOpenShellAlphaVar(int twoJ, HamiltonianMatrix* H)
{
    double alpha0 = Constant::Alpha;

    integrals->IncludeValenceSMS(false);
    std::string original_id = identifier;

    for(double x = -0.25; x <= 0.25; x += 0.125)
    {
        Constant::Alpha = alpha0 * sqrt(x+1.);
        Constant::AlphaSquared = alpha0 * alpha0 * (x+1.);

        *outstream << "\nx = " << x << std::endl;

        std::stringstream ss;
        ss << (int)(x * 8);
        identifier = original_id + '_' + ss.str();

        core->ToggleOpenShellCore();
        //core->Update();
        //excited->Update();
        //if(MBPT_CI)
        //    excited_mbpt->Update();
        //Write();
        Read();
        core->ToggleClosedShellCore();

        integrals->SetIdentifier(identifier);
        integrals->Update();

        H->GenerateMatrix();

        if(x == 0.)
            H->SolveMatrix(NumSolutions, twoJ, true);
        else
            H->SolveMatrix(NumSolutions, twoJ);
    }

    Constant::Alpha = alpha0;
    Constant::AlphaSquared = alpha0 * alpha0;

    identifier = original_id;
}
*/
