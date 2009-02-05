#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Atom.h"
#include "OutStreams.h"
#include "Universal/Constant.h"

#include "Basis/RStates.h"
#include "Basis/RSinStates.h"
#include "Basis/CustomBasis.h"
#include "Basis/ReadBasis.h"
#include "Basis/SubsetBasis.h"
#include "Basis/BSplineBasis.h"
#include "Basis/HartreeFockBasis.h"

#include "MBPT/CoreMBPTCalculator.h"
#include "Universal/ScalapackMatrix.h"
#include "Universal/ExpLattice.h"
#include "Universal/Eigensolver.h"
#include "RateCalculator.h"
#include "HartreeFock/ContinuumBuilder.h"
#include "HartreeFock/StateIntegrator.h"
#include <fstream>

#ifdef _MPI
    #ifdef _SCALAPACK
    #if !(_FUS)
        #define blacs_exit_ blacs_exit
    #endif
    extern "C"{
    void blacs_exit_(const int*);
    }
    #endif
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
        MPI::Intracomm& comm_world = MPI::COMM_WORLD; // Alias
        NumProcessors = comm_world.Get_size();
        ProcessorRank = comm_world.Get_rank();
    #else
        NumProcessors = 1;
        ProcessorRank = 0;
    #endif

    OutStreams::InitialiseStreams();

    try
    {   Atom A(90, 4, "ThIV", true);
        A.RunMultiple(true, true);
//        A.Run();
    }
    catch(std::bad_alloc& ba)
    {   *errstream << ba.what() << std::endl;
        exit(1);
    }

    #ifdef _MPI
        comm_world.Barrier();
        #ifdef _SCALAPACK
            // Continue should be non-zero, otherwise MPI_Finalise is called.
            int cont = 1;
            blacs_exit_(&cont);
        #endif

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
    //CreateHartreeFockBasis(false);
    CreateRBasis(true);
    //CreateBSplineBasis();

    mbpt = new CoreMBPTCalculator(lattice, core, excited_mbpt);
    //mbpt->UpdateIntegrals(excited);

    ConstStateIterator it = excited->GetConstStateIterator();
    while(!it.AtEnd())
    {
        StateInfo si = it.GetStateInfo();
        const DiscreteState* ds = it.GetState();

        if(excited->RetrieveSecondOrderSigma(si))
        {
            double dE = excited->GetSigmaMatrixElement(si);
            *outstream << (ds->Energy() + dE) * Constant::HartreeEnergy_cm << std::endl;

            DiscreteState brueckner = excited->GetStateWithSigma(si);
            *outstream << brueckner.Energy() * Constant::HartreeEnergy_cm << std::endl;
            *outstream << "\n" << brueckner.Name() << std::endl;
            for(unsigned int i = 0; i < 100; i+=10)
            {
                *outstream << brueckner.f[i]/ds->f[i] << "\t" << brueckner.g[i]/ds->g[i] << std::endl;
            }
        }

//        *outstream << ds->Name() << ": " << ds->Energy() * Constant::HartreeEnergy_cm
//                   << " + dE = " << (ds->Energy() + dE) * Constant::HartreeEnergy_cm << std::endl;
//        *outstream << (ds->Energy() + dE) * Constant::HartreeEnergy_cm << std::endl;
        it.Next();
    }
}

Atom::Atom(unsigned int atomic_number, int charge, const std::string& atom_identifier, bool read):
    Z(atomic_number), Charge(charge), identifier(atom_identifier),
    SD_CI(false), MBPT_CI(false), NumSolutions(6),
    excited(NULL), excited_mbpt(NULL), valence_mbpt(NULL),
    integrals(NULL), integralsMBPT(NULL), mbpt(NULL), sigma3(NULL)
{
    lattice = new Lattice(1000, 1.e-6, 50.);
    //        lattice->real_to_lattice(200.);
    //lattice = new ExpLattice(300, 1.e-5, 0.05);
    core = new Core(lattice, atomic_number, charge);
    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);
    //DebugOptions.HartreeEnergyUnits(true);

    if(read)
    {   Read();
        core->Update();
    }
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
    if(sigma3)
        delete sigma3;
    delete core;
    delete lattice;
}

void Atom::CreateBasis(bool UseMBPT)
{
    std::vector<unsigned int> num_states;
    num_states.push_back(1);
    num_states.push_back(1);
    num_states.push_back(1);
    num_states.push_back(1);
    num_states.push_back(1);

    if(UseMBPT)
    {   // Add extra states and waves to mbpt basis
        std::vector<unsigned int> num_states_mbpt;
        num_states_mbpt.push_back(2);
        num_states_mbpt.push_back(2);
        num_states_mbpt.push_back(3);
        num_states_mbpt.push_back(4);
        num_states_mbpt.push_back(4);
        num_states_mbpt.push_back(3);

        for(unsigned int i = 0; i < num_states_mbpt.size(); i++)
            num_states_mbpt[i] += 16;

        excited_mbpt->CreateExcitedStates(num_states_mbpt);        
    }

    excited->SetIdentifier(identifier);
    excited->CreateExcitedStates(num_states);
}

void Atom::CreateRBasis(bool UseMBPT)
{
    if(UseMBPT)
    {   excited_mbpt = new RSinStates(lattice, core);
        excited = new SubsetBasis(lattice, excited_mbpt);
    }
    else
    {   excited = new RSinStates(lattice, core);
    }

    CreateBasis(UseMBPT);
}

void Atom::CreateBSplineBasis(bool UseMBPT)
{
    BSplineBasis* basis;

    if(UseMBPT)
    {   excited_mbpt = new BSplineBasis(lattice, core);
        excited = new SubsetBasis(lattice, excited_mbpt);
        basis = dynamic_cast<BSplineBasis*>(excited_mbpt);
    }
    else
    {   excited = new BSplineBasis(lattice, core);
        basis = dynamic_cast<BSplineBasis*>(excited);
    }

    basis->SetParameters(40, 7, 50.);
    CreateBasis(UseMBPT);
}

void Atom::CreateBSplineBasis(double radius)
{
    excited = new BSplineBasis(lattice, core);
    dynamic_cast<BSplineBasis*>(excited)->SetParameters(40, 7, radius);

    CreateBasis(false);
}

void Atom::CreateCustomBasis(bool UseMBPT)
{
    if(UseMBPT)
    {   excited_mbpt = new CustomBasis(lattice, core);
        excited = new SubsetBasis(lattice, excited_mbpt);
    }
    else
    {   excited = new CustomBasis(lattice, core);
    }

    CreateBasis(UseMBPT);
}

void Atom::CreateHartreeFockBasis(bool UseMBPT)
{
    if(UseMBPT)
    {   excited_mbpt = new HartreeFockBasis(lattice, core);
        excited = new SubsetBasis(lattice, excited_mbpt);
    }
    else
    {   excited = new HartreeFockBasis(lattice, core);
    }

    CreateBasis(UseMBPT);
}

void Atom::Write() const
{
    if(ProcessorRank == 0)
    {
        std::string filename = identifier + ".core.atom";
        FILE* fp = fopen(filename.c_str(), "wb");

        // Output atomic data
        fwrite(&Z, sizeof(double), 1, fp);
        fwrite(&Charge, sizeof(double), 1, fp);

        // Output core states
        core->Write(fp);
        fclose(fp);

        if(excited)
        {   filename = identifier + ".excited.atom";
            fp = fopen(filename.c_str(), "wb");

            if(excited_mbpt)
                excited_mbpt->Write(fp);
            else
                excited->Write(fp);

            fclose(fp);
        }
    }
}

void Atom::Read()
{
    // Read core electron states
    std::string filename = identifier + ".core.atom";
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        exit(1);

    // Check that the stored ion is the same as this one!
    double stored_Z, stored_Charge;
    fread(&stored_Z, sizeof(double), 1, fp);
    fread(&stored_Charge, sizeof(double), 1, fp);
    if((stored_Z != Z) || (stored_Charge != Charge))
    {   fclose(fp);
        *errstream << "\nIncorrect stored state." << std::endl;
        exit(1);
    }

    core->Read(fp);
    fclose(fp);

    // Excited states
    if(excited)
    {   filename = identifier + ".excited.atom";
        fp = fopen(filename.c_str(), "rb");
        excited->Read(fp);

        if(excited_mbpt)
        {   fp = freopen(filename.c_str(), "rb", fp);
            excited_mbpt->Read(fp);
        }

        fclose(fp);
    }
}

void Atom::GenerateCowanInputFile()
{
    FILE* fp = fopen("test.txt", "wt");

    const DiscreteState* ds = excited->GetState(StateInfo(1, -1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(2, -1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(2, 1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(2, -2));
    PrintWavefunctionCowan(fp, ds);
/*
    ds = excited->GetState(StateInfo(4, 1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(4, -2));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(4, 3));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(4, -4));
    PrintWavefunctionCowan(fp, ds);
*/

    // Ground state
    Eigenstates* E = GetEigenstates(Symmetry(0, even));
    double energy_shift = -324.42474 - E->GetEigenvalues()[0];
//    double energy_shift = -1.7598437 - E->GetEigenvalues()[0];

    E = GetEigenstates(Symmetry(0, odd));
    E->PrintCowan(fp, energy_shift);

    E = GetEigenstates(Symmetry(2, odd));
    E->PrintCowan(fp, energy_shift);

    E = GetEigenstates(Symmetry(4, odd));
    E->PrintCowan(fp, energy_shift);

    fclose(fp);
}

void Atom::PrintWavefunctionCowan(FILE* fp, const DiscreteState* ds)
{
    // Extract A and B from expansion
    //   f(r) = A * r^g + B * r^(g + 1)
    // where g = |Kappa|.
    unsigned int r0 = 0;
    unsigned int r1 = 1;

    double R[4];
    R[0] = pow(lattice->R(r0), abs(ds->Kappa()));
    R[1] = pow(lattice->R(r0), abs(ds->Kappa()) + 1.);
    R[2] = pow(lattice->R(r1), abs(ds->Kappa()));
    R[3] = pow(lattice->R(r1), abs(ds->Kappa()) + 1.);

    double f[2];
    f[0] = ds->f[r0];
    f[1] = ds->f[r1];

    Eigensolver solver;
    if(!solver.SolveSimultaneousEquations(R, f, 2))
    {   *errstream << "PrintWavefunctionCowan: Can't get wavefunction expansion" << std::endl;
        exit(1);
    }

    // Upper component
    fprintf(fp, "     %4s     %5d%12.4E%12.4E%12.4E\n", ds->Name().c_str(), lattice->Size(), f[0], f[1], -ds->Energy());

    unsigned int count = 0;
    unsigned int i;
    for(i = 0; i < ds->Size(); i++)
    {   if(count == 5)
        {   fprintf(fp, "\n");
            count = 1;
        }
        else
            count++;
        fprintf(fp, "%14.7E", ds->f[i]);
    }
    while(i < lattice->Size())
    {   if(count == 5)
        {   fprintf(fp, "\n");
            count = 1;
        }
        else
            count++;
        fprintf(fp, "%14.7E", 0.0);
        i++;
    }
    
    // Lower component
    R[0] = pow(lattice->R(r0), abs(ds->Kappa()));
    R[1] = pow(lattice->R(r0), abs(ds->Kappa()) + 1.);
    R[2] = pow(lattice->R(r1), abs(ds->Kappa()));
    R[3] = pow(lattice->R(r1), abs(ds->Kappa()) + 1.);

    f[0] = ds->g[r0];
    f[1] = ds->g[r1];

    if(!solver.SolveSimultaneousEquations(R, f, 2))
    {   *errstream << "PrintWavefunctionCowan: Can't get wavefunction expansion" << std::endl;
        exit(1);
    }

    fprintf(fp, "\n     %4s     %5d%12.4E%12.4E\n", ds->Name().c_str(), lattice->Size(), f[0], f[1]);

    count = 0;
    for(i = 0; i < ds->Size(); i++)
    {   if(count == 5)
        {   fprintf(fp, "\n");
            count = 1;
        }
        else
            count++;
        fprintf(fp, "%14.7E", ds->g[i]*Constant::Alpha);
    }
    while(i < lattice->Size())
    {   if(count == 5)
        {   fprintf(fp, "\n");
            count = 1;
        }
        else
            count++;
        fprintf(fp, "%14.7E", 0.0);
        i++;
    }
    
    fprintf(fp, "\n");
}
