#include "Include.h"
#include "Atom.h"
#include "OutStreams.h"
#include "Basis/RStates.h"
#include "Basis/RSinStates.h"
#include "Basis/CustomBasis.h"
#include "Universal/Constant.h"
#include "Basis/BSplineBasis.h"
#include "MBPT/MBPTCalculator.h"

#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"

#ifdef _MPI
#include <mpi.h>
#endif
// MPI details (if not used, we can have NumProcessors == 1)
unsigned int NumProcessors;
unsigned int ProcessorRank;

// The debug options for the whole program.
Debug DebugOptions;

int main(int argc, char* argv[])
{
    #ifdef _MPI
        MPI::Init(argc, argv);
        MPI::Intracomm comm_world = MPI::COMM_WORLD; // Alias
        NumProcessors = comm_world.Get_size();
        ProcessorRank = comm_world.Get_rank();
    #else
        NumProcessors = 1;
        ProcessorRank = 0;
    #endif

    OutStreams::InitialiseStreams();

    try
    {   Atom A(12, 2, "MgI001");
        A.RunOpen();
    }
    catch(std::bad_alloc& ba)
    {   *errstream << ba.what() << std::endl;
        exit(1);
    }

    #ifdef _MPI
        MPI::Finalize();
    #endif

    *outstream << "\nFinished" << std::endl;
    OutStreams::FinaliseStreams();

    PAUSE
    return 0;
}

void Atom::Run()
{
    DebugOptions.LogFirstBuild(false);
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);

    //CreateCustomBasis();
    //CreateRBasis();
    CreateBSplineBasis();

    DebugOptions.OutputHFExcited(false);

    DoClosedShellSMS(true);
}

Atom::Atom(unsigned int atomic_number, int charge, const std::string& atom_identifier, bool read):
    Z(atomic_number), Charge(charge), identifier(atom_identifier),
    SD_CI(false), MBPT_CI(false), NumSolutions(6), excited(NULL), excited_mbpt(NULL),
    integrals(NULL), integralsMBPT(NULL), mbpt(NULL)
{
    lattice = new Lattice(1000, 1.e-6, 50.);
    core = new Core(lattice, atomic_number, charge);
    //DebugOptions.LogFirstBuild(true);
    //DebugOptions.LogHFIterations(true);
    //DebugOptions.HartreeEnergyUnits(true);
    if(read)
        Read();
    else
        core->Initialise();
}

Atom::~Atom(void)
{
    if(integrals)
        delete integrals;
    if(excited)
        delete excited;
    if(excited_mbpt)
        delete excited_mbpt;
    if(mbpt)
        delete mbpt;
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
        *errstream << "\nIncorrect stored state." << std::endl;
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
    DiscreteState ds = excited->GetStateWithSigma(info);
    return ds.Energy();
}

void Atom::CreateRBasis(const StateInfo* ionised)
{
    excited = new RSinStates(lattice, core);
    excited->SetIdentifier(&identifier);

    if(ionised)
        core->Ionise(*ionised);

    std::vector<unsigned int> num_states;
    num_states.push_back(20);
    num_states.push_back(20);
    num_states.push_back(20);
    num_states.push_back(20);
    num_states.push_back(20);

    excited->CreateExcitedStates(num_states);
}

void Atom::CreateBSplineBasis(const StateInfo* ionised)
{
    excited = new BSplineBasis(lattice, core);
    excited->SetIdentifier(&identifier);

    dynamic_cast<BSplineBasis*>(excited)->SetParameters(40, 7, 45.);

    if(ionised)
        core->Ionise(*ionised);

    std::vector<unsigned int> num_states;
    num_states.push_back(7);
    num_states.push_back(7);
    num_states.push_back(7);
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

void Atom::DoClosedShellSMS(bool include_mbpt)
{
    MBPTCalculator mbpt(lattice, core, excited);
    const unsigned int max_k = 3;
    double totals[max_k];

    for(double ais = -0.002; ais <= 0.002; ais += 0.001)
    {
        core->SetNuclearInverseMass(ais);
        core->Update();
        excited->Update();

        const State* ds;
        int k2;
        for(k2 = 1; k2 <= max_k; k2++)
        {
            int kappa = k2/2;
            if(k2%2)
                kappa = -kappa-1;

            if(k2 <= 3)
                ds = excited->GetState(StateInfo(2, kappa));
            else
                ds = excited->GetState(StateInfo(5, kappa));

            if(include_mbpt)
            {   totals[k2-1] = mbpt.GetSecondOrderSigma(ds);
            }
            else
                totals[k2-1] = ds->Energy();
        }

        *outstream << "\nNuclearInverseMass = " << ais << std::endl;
        for(k2 = 0; k2 < max_k; k2++)
            *outstream << std::setprecision(15) << totals[k2]*Constant::HartreeEnergy_cm << std::endl;
        *outstream << std::endl;
    }
}

void Atom::DoClosedShellVolumeShift(bool include_mbpt)
{
    core->CalculateVolumeShiftPotential(0.05/Constant::AtomicToFermi);

    MBPTCalculator mbpt(lattice, core, excited);
    const unsigned int max_k = 5;
    double totals[max_k];

    for(double ais = -100.; ais <= 100.; ais += 50.)
    {
        core->SetVolumeShiftParameter(ais);
        core->Update();
        excited->Update();

        const State* ds;
        int k2;
        for(k2 = 1; k2 <= max_k; k2++)
        {
            int kappa = k2/2;
            if(k2%2)
                kappa = -kappa-1;

            if(k2 <= 3)
                ds = excited->GetState(StateInfo(5, kappa));
            else
                ds = excited->GetState(StateInfo(5, kappa));

            if(include_mbpt)
                totals[k2-1] = mbpt.GetSecondOrderSigma(ds);
            else
                totals[k2-1] = ds->Energy();
        }

        *outstream << "\nVolumeShiftParameter = " << ais << std::endl;
        for(k2 = 0; k2 < max_k; k2++)
            *outstream << std::setprecision(15) << totals[k2]*Constant::HartreeEnergy_cm << std::endl;
        *outstream << std::endl;
    }

    core->SetVolumeShiftParameter(0.);
}

void Atom::DoClosedShellAlphaVar(bool include_mbpt)
{
    MBPTCalculator mbpt(lattice, core, excited);
    const unsigned int max_k = 5;
    double totals[max_k];

    double alpha0 = Constant::Alpha;

    for(double x = -0.25; x <= 0.25; x += 0.125)
    {
        Constant::Alpha = alpha0 * sqrt(x+1.);
        Constant::AlphaSquared = alpha0 * alpha0 * (x+1.);

        core->Update();
        excited->Update();

        const State* ds;
        int k2;
        for(k2 = 1; k2 <= max_k; k2++)
        {
            int kappa = k2/2;
            if(k2%2)
                kappa = -kappa-1;

            if(k2 <= 3)
                ds = excited->GetState(StateInfo(5, kappa));
            else
                ds = excited->GetState(StateInfo(5, kappa));

            if(include_mbpt)
                totals[k2-1] = mbpt.GetSecondOrderSigma(ds);
            else
                totals[k2-1] = ds->Energy();
        }

        *outstream << "\nx = " << x << std::endl;
        for(k2 = 0; k2 < max_k; k2++)
            *outstream << totals[k2]*Constant::HartreeEnergy_cm << std::endl;
        *outstream << std::endl;
    }

    Constant::Alpha = alpha0;
    Constant::AlphaSquared = alpha0 * alpha0;
}
