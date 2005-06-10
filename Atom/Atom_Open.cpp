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
    MBPT_CI = true;

    if(MBPT_CI)
    {   integrals = new CIIntegralsMBPT(*excited);
        integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);

        excited_mbpt = new BSplineBasis(lattice, core);
        excited_mbpt->SetIdentifier(&identifier);
        dynamic_cast<BSplineBasis*>(excited_mbpt)->SetParameters(40, 7, 45.);
        std::vector<unsigned int> num_states_per_l;
        num_states_per_l.push_back(30);
        num_states_per_l.push_back(30);
        num_states_per_l.push_back(30);
        num_states_per_l.push_back(29);
        num_states_per_l.push_back(28);
        excited_mbpt->CreateExcitedStates(num_states_per_l);

        mbpt = new MBPTCalculator(lattice, core, excited_mbpt);

        integralsMBPT->IncludeMBPT1(true, mbpt);
        integralsMBPT->IncludeMBPT2(true, mbpt);
    }
    else
    {   integrals = new CIIntegrals(*excited);
    }

    Configuration config;
    config.SetOccupancy(NonRelInfo(3, 0), 2);

    unsigned int two_j;
    NumSolutions = 6;

    *outstream << "\nGS Parity:\n" << std::endl;
    for(two_j = 0; two_j <= 6; two_j += 2)
    {
        RelativisticConfigList rlist;
        HamiltonianMatrix* H = CreateHamiltonian(two_j, config, rlist);

        //OpenShellEnergy(two_j, H);
        DoOpenShellSMS(two_j, H);
        //DoOpenShellVolumeShift(two_j, H);

        delete H;
    }

    config.RemoveSingleParticle(NonRelInfo(3, 0));
    config.AddSingleParticle(NonRelInfo(3, 1));

    NumSolutions = 6;

    *outstream << "\nOpposite Parity:\n" << std::endl;
    for(two_j = 0; two_j <= 4; two_j+=2)
    {
        RelativisticConfigList rlist;
        HamiltonianMatrix* H = CreateHamiltonian(two_j, config, rlist);
        
        //OpenShellEnergy(two_j, H);
        DoOpenShellSMS(two_j, H);
        //DoOpenShellVolumeShift(two_j, H);

        delete H;
    }
}

HamiltonianMatrix* Atom::CreateHamiltonian(int twoJ, const Configuration& config, RelativisticConfigList& rlist)
{
    unsigned int electron_excitations;
    if(SD_CI)
        electron_excitations = 2;
    else
        electron_excitations = config.NumParticles();

    ConfigList nrlist;
    nrlist.push_back(config);

    rlist.clear();

    ConfigGenerator generator(excited);
    generator.GenerateMultipleExcitations(nrlist, electron_excitations);
    generator.GenerateRelativisticConfigs(nrlist, rlist);
    generator.GenerateProjections(rlist, twoJ);

    HamiltonianMatrix* H;

    #ifdef _MPI
        H = new MPIHamiltonianMatrix(*integrals, rlist);
    #else
        H = new HamiltonianMatrix(*integrals, rlist);
    #endif

    return H;
}

void Atom::CheckMatrixSizes()
{
    // Two electron integral storage size
    *outstream << "Num coulomb integrals: " << integrals->GetStorageSize() << std::endl;

    SD_CI = true;
    Configuration config;
    config.SetOccupancy(NonRelInfo(2, 0), 2);

    unsigned int two_j;

    *outstream << "\nGS Parity:\n" << std::endl;
    for(two_j = 0; two_j <= 0; two_j += 6)
    {
        *outstream << "J = " << two_j/2. << std::endl;
        OpenShellEnergy(two_j, config, true);
    }

    config.RemoveSingleParticle(NonRelInfo(2, 0));
    config.AddSingleParticle(NonRelInfo(2, 1));

    *outstream << "\nOpposite Parity:\n" << std::endl;
    for(two_j = 0; two_j <= 4; two_j+=2)
    {
        *outstream << "J = " << two_j/2. << std::endl;
        OpenShellEnergy(two_j, config, true);
    }
}

void Atom::OpenShellEnergy(int twoJ, const Configuration& config, bool size_only)
{
    if(size_only)
    {   // We really don't want to allocate space for the HamiltonianMatrix.
        unsigned int electron_excitations;
        if(SD_CI)
            electron_excitations = 2;
        else
            electron_excitations = config.NumParticles();

        ConfigList nrlist;
        nrlist.push_back(config);

        RelativisticConfigList rlist;
        ConfigGenerator generator(excited);
        generator.GenerateMultipleExcitations(nrlist, electron_excitations);
        generator.GenerateRelativisticConfigs(nrlist, rlist);
        *outstream << " Number of non-rel configurations = " << nrlist.size() << std::endl;
        *outstream << " Number of rel configurations = " << rlist.size() << std::endl;
        generator.GenerateProjections(rlist, twoJ);
        
        unsigned int N = 0;
        unsigned int i = 0;
        RelativisticConfigList::const_iterator it = rlist.begin();
        while(it != rlist.end())
        {   N += it->NumJStates();
            it++;
        }
        *outstream << " Number of J-configurations = " << N << std::endl;
    }
    else
    {   RelativisticConfigList rlist;
        HamiltonianMatrix* H = CreateHamiltonian(twoJ, config, rlist);

        OpenShellEnergy(twoJ, H);

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
}

void Atom::DoOpenShellSMS(int twoJ, HamiltonianMatrix* H)
{
    integrals->IncludeValenceSMS(true);

    std::string original_id = identifier;
    
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
        //Write();
        Read();
        if(MBPT_CI)
            excited_mbpt->Update();
        core->ToggleClosedShellCore();

        if(NumProcessors > 1)
        {   std::stringstream proc;
            proc << ProcessorRank;
            std::string id = identifier + '_' + proc.str();
            integrals->SetIdentifier(id);
        }
        else
            integrals->SetIdentifier(identifier);

        integrals->Update();

        //integralsMBPT->WriteSigmaPotentials();
        //integrals->WriteOneElectronIntegrals();
        //integrals->WriteTwoElectronIntegrals();

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

void Atom::DoOpenShellVolumeShift(int twoJ, HamiltonianMatrix* H)
{
    static bool calculated_potential = false;
    if(!calculated_potential)
    {   core->CalculateVolumeShiftPotential(0.05/Constant::AtomicToFermi);
        calculated_potential = true;
    }

    integrals->IncludeValenceSMS(false);
    for(double ais = -100.; ais <= 100.; ais += 50.)
    {
        *outstream << "\nVolumeShiftParameter = " << std::setprecision(3) << ais << std::endl;

        core->ToggleOpenShellCore();
        core->SetVolumeShiftParameter(ais);
        core->Update();
        excited->Update();

        core->ToggleClosedShellCore();
        integrals->Update();
        H->GenerateMatrix();

        if(ais == 0.0)  // Get g-factors
            H->SolveMatrix(NumSolutions, twoJ, true);
        else
            H->SolveMatrix(NumSolutions, twoJ);
    }
}

void Atom::DoOpenShellAlphaVar(int twoJ, HamiltonianMatrix* H)
{
    double alpha0 = Constant::Alpha;
    integrals->IncludeValenceSMS(false);

    for(double x = -0.25; x <= 0.25; x += 0.125)
    {
        Constant::Alpha = alpha0 * sqrt(x+1.);
        Constant::AlphaSquared = alpha0 * alpha0 * (x+1.);

        *outstream << "\nx = " << x << std::endl;

        core->ToggleOpenShellCore();
        core->Update();
        excited->Update();

        core->ToggleClosedShellCore();
        integrals->Update();
        H->GenerateMatrix();

        if(x == 0.)
            H->SolveMatrix(NumSolutions, twoJ, true);
        else
            H->SolveMatrix(NumSolutions, twoJ);
    }

    Constant::Alpha = alpha0;
    Constant::AlphaSquared = alpha0 * alpha0;
}
