#include "Include.h"
#include "Atom.h"
#include "Basis/HFExcitedStates.h"
#include "Basis/RStates.h"
#include "Basis/RSinStates.h"
#include "Basis/CustomBasis.h"
#include "Universal/Constant.h"
#include "Basis/BSplineBasis.h"
#include "MBPT/MBPTCalculator.h"
#include <fstream>

int main(int argc, char* argv[])
{
#ifdef UNIX   // Don't output cerr to screen
    std::ofstream error_file("error.txt");
    std::cerr.rdbuf(error_file.rdbuf());
#endif

    try
    {   Atom A(26, 2, "FeII");
        A.RunOpen();
    }
    catch(std::bad_alloc& ba)
    {   std::cout << ba.what() << std::endl;
        exit(1);
    }

#ifndef UNIX
#ifdef _DEBUG
    std::cout << "\nFinished" << std::endl;
    getchar();
#endif
#endif

    return 0;
}

void Atom::Run()
{
    GetDebugOptions().DebugFirstBuild(false);
    GetDebugOptions().DebugHFIterations(false);
    GetDebugOptions().DebugHFExcited(true);
    GetDebugOptions().HartreeEnergyUnits(true);

    //CreateHFBasis();
    //CreateCustomBasis();
    //CreateRBasis();
    CreateBSplineBasis();

    GetDebugOptions().DebugHFExcited(false);

    DoClosedShellAlphaVar();
}

Atom::Atom(unsigned int atomic_number, int charge, const std::string& atom_identifier, bool read):
    Z(atomic_number), Charge(charge), identifier(atom_identifier),
    SD_CI(false), NumSolutions(6), excited(NULL)
{
    lattice = new Lattice(1000, 1.e-6, 50.);
    core = new Core(lattice, atomic_number, charge);
    //GetDebugOptions().DebugFirstBuild(true);
    //GetDebugOptions().DebugHFIterations(true);
    //GetDebugOptions().HartreeEnergyUnits(true);
    if(read)
        Read();
    else
        core->Initialise();
}

Atom::~Atom(void)
{
    if(excited)
        delete excited;
    delete core;
    delete lattice;
}

void Atom::Write() const
{
    std::string filename = identifier + ".atom";
    FILE* fp = fopen(filename.c_str(), "wb");

    // Output atomic data
    fwrite(&Z, sizeof(double), 1, fp);
    fwrite(&Charge, sizeof(double), 1, fp);

    // Output electron states
    core->Write(fp);
    if(excited)
        excited->Write(fp);
    else
    {   unsigned int zero = 0;
        fwrite(&zero, sizeof(unsigned int), 1, fp);
    }
    fclose(fp);
}

void Atom::Read()
{
    std::string filename = identifier + ".atom";
    FILE* fp = fopen(filename.c_str(), "rb");

    // Check that the stored ion is the same as this one!
    double stored_Z, stored_Charge;
    fread(&stored_Z, sizeof(double), 1, fp);
    fread(&stored_Charge, sizeof(double), 1, fp);
    if((stored_Z != Z) || (stored_Charge != Charge))
    {   fclose(fp);
        std::cerr << "\nIncorrect stored state." << std::endl;
        exit(1);
    }

    // Read electron states
    core->Read(fp);
    if(excited)
        excited->Read(fp);

    fclose(fp);
}

double Atom::GetEnergy(const StateInfo& info)
{
    if(!info.Discrete())
        return 0.;

    DiscreteState ds = excited->GetStateWithSigma(info);
    return ds.Energy();
}

void Atom::CreateHFBasis()
{
    excited = new HFExcitedStates(lattice, core);
    excited->SetIdentifier(&identifier);

    std::vector<unsigned int> num_states;
    num_states.push_back(1);
    num_states.push_back(1);
    num_states.push_back(1);

    excited->CreateExcitedStates(num_states);
    // dynamic_cast<HFExcitedStates*>(excited)->CreateContinuum(0.1, 4.2, 25, 7);
    // CreateContinuum(0.1, 4.2, 25, 11)
}

void Atom::CreateRBasis(const StateInfo* ionised)
{
    excited = new RSinStates(lattice, core);
    excited->SetIdentifier(&identifier);

    if(ionised)
        core->Ionise(*ionised);

    std::vector<unsigned int> num_states;
    num_states.push_back(2);
    num_states.push_back(2);
    num_states.push_back(2);
    //num_states.push_back(13);
    //num_states.push_back(13);
    //num_states.push_back(12);
    //num_states.push_back(7);
    
    //{3, 3, 4, 4};{5, 5, 6, 2};{6, 6, 3};
    excited->CreateExcitedStates(num_states);
}

void Atom::CreateBSplineBasis(const StateInfo* ionised)
{
    excited = new BSplineBasis(lattice, core);
    excited->SetIdentifier(&identifier);

    dynamic_cast<BSplineBasis*>(excited)->SetParameters(40, 7, 50.);

    if(ionised)
        core->Ionise(*ionised);

    std::vector<unsigned int> num_states;
    num_states.push_back(10);
    num_states.push_back(10);
    num_states.push_back(10);
    num_states.push_back(6);

    excited->CreateExcitedStates(num_states);
}

void Atom::CreateCustomBasis(const StateInfo* ionised)
{
    excited = new CustomBasis(lattice, core);
    excited->SetIdentifier(&identifier);

    if(ionised)
        core->Ionise(*ionised);

    std::vector<unsigned int> num_states;
    excited->CreateExcitedStates(num_states);
}

void Atom::GetSigma(const StateInfo& info)
{
    DiscreteState* ds = dynamic_cast<DiscreteState*>(excited->GetState(info));
    if(ds)
    {   MBPTCalculator MC(lattice, core, excited);
        printf("%s %.5f\n\n", ds->Name().c_str(), MC.GetSecondOrderSigma(ds)*Constant::HartreeEnergy_cm);
    }
}

void Atom::DoClosedShellSMS(bool include_mbpt)
{
    MBPTCalculator mbpt(lattice, core, excited);
    double totals[5];

    for(double ais = -0.002; ais <= 0.002; ais += 0.001)
    {
        core->SetNuclearInverseMass(ais);
        core->Update();
        excited->Update();

        const State* ds;
        int k2;
        for(k2 = 1; k2 <= 5; k2++)
        {
            int kappa = k2/2;
            if(k2%2)
                kappa = -kappa-1;

            if(k2 <=3)
                ds = excited->GetState(StateInfo(6, kappa));
            else
                ds = excited->GetState(StateInfo(5, kappa));

            if(include_mbpt)
                totals[k2-1] = mbpt.GetSecondOrderSigma(ds);
            else
                totals[k2-1] = ds->Energy();
        }

        std::cout << "\nNuclearInverseMass = " << ais << std::endl;
        for(k2 = 0; k2 < 5; k2++)
            printf("%.5f\n", totals[k2]*Constant::HartreeEnergy_cm);
        std::cout << std::endl;
    }
}

void Atom::DoClosedShellVolumeShift(bool include_mbpt)
{
    core->CalculateVolumeShiftPotential(0.01/Constant::AtomicToFermi);

    MBPTCalculator mbpt(lattice, core, excited);
    double totals[5];

    for(double ais = -100.; ais <= 100.; ais += 20.)
    {
        core->SetVolumeShiftParameter(ais);
        core->Update();
        excited->Update();

        const State* ds;
        int k2;
        for(k2 = 1; k2 <= 5; k2++)
        {
            int kappa = k2/2;
            if(k2%2)
                kappa = -kappa-1;

            if(k2 <= 3)
                ds = excited->GetState(StateInfo(6, kappa));
            else
                ds = excited->GetState(StateInfo(5, kappa));

            if(include_mbpt)
                totals[k2-1] = mbpt.GetSecondOrderSigma(ds);
            else
                totals[k2-1] = ds->Energy();
        }

        std::cout << "\nVolumeShiftParameter = " << ais << std::endl;
        for(k2 = 0; k2 < 5; k2++)
            printf("%.5f\n", totals[k2]*Constant::HartreeEnergy_cm);
        std::cout << std::endl;
    }

    core->SetVolumeShiftParameter(0.);
}

void Atom::DoClosedShellAlphaVar(bool include_mbpt)
{
    MBPTCalculator mbpt(lattice, core, excited);
    double totals[5];

    const unsigned int N = 3;

    double alpha0 = Constant::Alpha;

    for(double x = -0.25; x <= 0.25; x += 0.125)
    {
        Constant::Alpha = alpha0 * sqrt(x+1.);
        Constant::AlphaSquared = alpha0 * alpha0 * (x+1.);

        core->Update();
        excited->Update();

        const State* ds;
        int k2;
        for(k2 = 1; k2 <= N; k2++)
        {
            int kappa = k2/2;
            if(k2%2)
                kappa = -kappa-1;

            if(k2 <= 3)
                ds = excited->GetState(StateInfo(2, kappa));
            else
                ds = excited->GetState(StateInfo(3, kappa));

            if(include_mbpt)
                totals[k2-1] = mbpt.GetSecondOrderSigma(ds);
            else
                totals[k2-1] = ds->Energy();
        }

        std::cout << "\nx = " << x << std::endl;
        for(k2 = 0; k2 < N; k2++)
            printf("%.5f\n", totals[k2]*Constant::HartreeEnergy_cm);
        std::cout << std::endl;
    }

    Constant::Alpha = alpha0;
    Constant::AlphaSquared = alpha0 * alpha0;
}

