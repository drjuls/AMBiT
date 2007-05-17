#include "Include.h"
#include "Atom.h"
#include "OutStreams.h"
#include "Basis/RStates.h"
#include "Basis/RSinStates.h"
#include "Basis/CustomBasis.h"
#include "Universal/Constant.h"
#include "Basis/BSplineBasis.h"
#include "MBPT/MBPTCalculator.h"
#include "Universal/ScalapackMatrix.h"

#include "Universal/ExpLattice.h"
#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"
#include "Universal/Eigensolver.h"

#include "RateCalculator.h"

#ifdef _MPI
    #ifdef _SCALAPACK
    #if !(_FUS)
        #define blacs_exit_ blacs_exit
    #endif
    extern "C"{
    void blacs_exit_(const int*);
    }
    #endif
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
        MPI::Intracomm& comm_world = MPI::COMM_WORLD; // Alias
        NumProcessors = comm_world.Get_size();
        ProcessorRank = comm_world.Get_rank();
    #else
        NumProcessors = 1;
        ProcessorRank = 0;
    #endif

    OutStreams::InitialiseStreams();

    try
    {   Atom A(6, 4, "CIII");
        A.RunOpen();
//        A.GenerateCowanInputFile();
        RateCalculator rate(A.GetBasis());
        double A1 = rate.CalculateAugerRate(&A, Symmetry(2, even), &A, Symmetry(2, odd), 4);
        *outstream << A1 << std::endl;
//        rate.CalculateAllDipoleStrengths(&A, Symmetry(2, odd), 2);
//        rate.CalculateDipoleStrength(&A, Symmetry(0, even), 0, Symmetry(2, odd), 2);
//        rate.CalculateDipoleStrength(&A, Symmetry(2, odd), 0, Symmetry(2, even), 1);
    }
    catch(std::bad_alloc& ba)
    {   *errstream << ba.what() << std::endl;
        exit(1);
    }

    #ifdef _MPI
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
    //CreateRBasis();
    //CreateBSplineBasis();

    DebugOptions.OutputHFExcited(false);

    ContinuumState* cs = new ContinuumState(lattice, 1.0, -1, core->GetHFPotential().size());
    unsigned int loops = core->CalculateExcitedState(cs);
    *outstream << loops << "\n" << std::endl;

    *outstream << cs->Size() << std::endl;
    for(unsigned int i=0; i<cs->Size(); i+=1)
        *outstream << i << "\t " << cs->f[i] << std::endl;
    //DoClosedShellSMS(true);
}

Atom::Atom(unsigned int atomic_number, int charge, const std::string& atom_identifier, bool read):
    Z(atomic_number), Charge(charge), identifier(atom_identifier),
    SD_CI(false), MBPT_CI(false), NumSolutions(6),
    excited(NULL), excited_mbpt(NULL),
    integrals(NULL), integralsMBPT(NULL), mbpt(NULL), sigma3(NULL)
{
    lattice = new Lattice(1500, 1.e-6, 200.);
    //lattice = new Lattice(1000, 1.e-6, 50.);
    //lattice = new ExpLattice(350, 1.e-5, 0.05);
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
    if(sigma3)
        delete sigma3;
    delete core;
    delete lattice;
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
    num_states.push_back(6);
    num_states.push_back(6);
    num_states.push_back(7);
    num_states.push_back(6);

    excited->CreateExcitedStates(num_states);
}

void Atom::CreateBSplineBasis(const StateInfo* ionised)
{
    excited = new BSplineBasis(lattice, core);
    excited->SetIdentifier(&identifier);

    dynamic_cast<BSplineBasis*>(excited)->SetParameters(40, 7, 25.);

    if(ionised)
        core->Ionise(*ionised);

    std::vector<unsigned int> num_states;
    num_states.push_back(7);
    num_states.push_back(7);
    num_states.push_back(6);
    num_states.push_back(5);
//    num_states.push_back(13);

    excited->CreateExcitedStates(num_states);
}

void Atom::CreateCustomBasis(const StateInfo* ionised)
{
    excited = new CustomBasis(lattice, core);
    excited->SetIdentifier(&identifier);

    if(ionised)
        core->Ionise(*ionised);

    std::vector<unsigned int> num_states;
    num_states.push_back(3);
    num_states.push_back(3);
    num_states.push_back(4);
    num_states.push_back(3);

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

void Atom::DoClosedShellFSModifyR(bool include_mbpt)
{
    MBPTCalculator mbpt(lattice, core, excited);
    const unsigned int max_k = 5;
    double totals[max_k];

    *outstream << "\nThickness = " << std::setprecision(3) << core->GetNuclearThickness()*Constant::AtomicToFermi;

    for(double r = 0.; r <= 12.; r += 3.)
    {
        core->SetNuclearRadius(r/Constant::AtomicToFermi);
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
                ds = excited->GetState(StateInfo(4, kappa));
            else
                ds = excited->GetState(StateInfo(3, kappa));

            if(include_mbpt)
                totals[k2-1] = mbpt.GetSecondOrderSigma(ds);
            else
                totals[k2-1] = ds->Energy();
        }

        *outstream << "\nRadius = " << r << std::endl;
        for(k2 = 0; k2 < max_k; k2++)
            *outstream << std::setprecision(15) << totals[k2]*Constant::HartreeEnergy_cm << std::endl;
        *outstream << std::endl;
    }
}

void Atom::DoClosedShellVolumeShift(bool include_mbpt)
{
    double deltaR = 0.05;
    core->CalculateVolumeShiftPotential(deltaR/Constant::AtomicToFermi);

    MBPTCalculator mbpt(lattice, core, excited);
    const unsigned int max_k = 5;
    double totals[max_k];

    *outstream << "\nThickness = " << std::setprecision(3) << core->GetNuclearThickness()*Constant::AtomicToFermi;
    *outstream << "\nRadius    = " << std::setprecision(3) << core->GetNuclearRadius()*Constant::AtomicToFermi;
    *outstream << "\nd(Radius) = " << std::setprecision(3) << deltaR;

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

void Atom::GenerateCowanInputFile()
{
    FILE* fp = fopen("test.txt", "wt");

    const DiscreteState* ds = core->GetState(StateInfo(1, -1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(2, -1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(2, 1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(2, -2));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(4, 1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(4, -2));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(4, 3));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(4, -4));
    PrintWavefunctionCowan(fp, ds);

    Eigenstates* E = GetEigenstates(Symmetry(0, even));
    double energy_shift = -1.7598437 - E->GetEigenvalues()[0];
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
