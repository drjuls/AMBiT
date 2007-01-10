#include "Include.h"
#include "Atom.h"
#include "Universal/Constant.h"
#include "HartreeFock/NonRelInfo.h"
#include "Configuration/ConfigGenerator.h"
#include "Configuration/HamiltonianMatrix.h"
#include "Configuration/MPIHamiltonianMatrix.h"
#include "Basis/BSplineBasis.h"

void Atom::RunOpen()
{
    DebugOptions.LogFirstBuild(false);
    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);

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
        dynamic_cast<BSplineBasis*>(excited_mbpt)->SetParameters(40, 7, 45.);
        std::vector<unsigned int> num_states_per_l;
        num_states_per_l.push_back(30);
        num_states_per_l.push_back(30);
        num_states_per_l.push_back(31);
        num_states_per_l.push_back(30);
        num_states_per_l.push_back(29);
        excited_mbpt->CreateExcitedStates(num_states_per_l);

        mbpt = new MBPTCalculator(lattice, core, excited_mbpt);

        integralsMBPT->IncludeMBPT1(true, mbpt);
        integralsMBPT->IncludeMBPT2(true, mbpt);
        integralsMBPT->IncludeExtraBoxDiagrams(true);
        integralsMBPT->SetExtraBoxDiagramLimits(4, 4);
    }
    else
    {   integrals = new CIIntegrals(*excited);
        integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);
    }

    //integralsMBPT->SetTwoElectronStorageLimits(4, 4);
    //sigma3 = new Sigma3Calculator(lattice, core, excited);

    GenerateFromFile = true;
    generator = new ConfigFileGenerator(excited);

    unsigned int NumParticles = 3;
    //CheckMatrixSizes();

    unsigned int two_j;
    NumSolutions = 3;

    *outstream << "\nGS Parity:\n" << std::endl;
    for(two_j = 3; two_j <= 9; two_j += 6)
    {
        HamiltonianMatrix* H = CreateHamiltonian(two_j, NumParticles, even);
        //H->IncludeSigma3(sigma3);

        //OpenShellEnergy(two_j, H);
        DoOpenShellSMS(two_j, H);
        //DoOpenShellVolumeShift(two_j, H);

        delete H;
    }
    
    NumSolutions = 4;

    *outstream << "\nOpposite Parity:\n" << std::endl;
    for(two_j = 9; two_j <= 11; two_j+=2)
    {
        HamiltonianMatrix* H = CreateHamiltonian(two_j, NumParticles, odd);
        //H->IncludeSigma3(sigma3);
        
        //OpenShellEnergy(two_j, H);
        DoOpenShellSMS(two_j, H);
        //DoOpenShellVolumeShift(two_j, H);

        delete H;
    }
}

HamiltonianMatrix* Atom::CreateHamiltonian(int twoJ, unsigned int num_particles, Parity P)
{
    unsigned int electron_excitations;
    if(SD_CI)
        electron_excitations = 2;
    else
        electron_excitations = num_particles;

    std::set<Configuration> leading_configs;

    Configuration config1;
    config1.SetOccupancy(NonRelInfo(3, 2), 2);
    config1.SetOccupancy(NonRelInfo(4, 0), 1);
    leading_configs.insert(config1);

    Configuration config2;
    config2.SetOccupancy(NonRelInfo(3, 2), 2);
    config2.SetOccupancy(NonRelInfo(4, 1), 1);
    leading_configs.insert(config2);

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

    unsigned int NumParticles = 3;
    unsigned int two_j;

    *outstream << "\nGS Parity:\n" << std::endl;
    for(two_j = 3; two_j <= 9; two_j += 6)
    {
        *outstream << "J = " << two_j/2. << std::endl;
        HamiltonianMatrix* H = CreateHamiltonian(two_j, NumParticles, even);
        *outstream << " Number of rel configurations = " << generator->GetRelConfigs()->size() << std::endl;
        delete H;
    }

    *outstream << "\nOpposite Parity:\n" << std::endl;
    for(two_j = 9; two_j <= 11; two_j+=2)
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
    core->SetVolumeShiftParameter(0.);
    core->Update();
    excited->Update();
    core->ToggleClosedShellCore();
    
    integrals->IncludeValenceSMS(false);
    integrals->Update();
    H->GenerateMatrix();
    //H->PollMatrix();
    H->SolveMatrix(NumSolutions, twoJ, true);

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
