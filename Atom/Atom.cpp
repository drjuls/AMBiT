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

Atom::Atom(GetPot userInput, unsigned int atomic_number, int num_electrons, const std::string& atom_identifier):
    userInput_(userInput), Z(atomic_number), identifier(atom_identifier)
{
    Charge = Z - num_electrons;

    lattice = NULL;
    core = NULL;
    excited = NULL;
    excited_mbpt = NULL;
    integrals = NULL;
    integralsMBPT = NULL;
    mbpt = NULL;
    valence_mbpt = NULL;
    sigma3 = NULL;

    // Operational parameters
    useRead = true;
    useWrite = true;

    // CI + MBPT parameters
    NumSolutions = 6;
    MaxEnergy = 0.0;
    mbptBasisString = "";
    ciBasisString = "";
    check_size_only = false;
    save_eigenstates = false;
    generate_mbpt_integrals = false;
    includeSigma1 = false;
    includeSigma2 = false;
    includeSigma3 = false;

    // Multiple run parameters
    multiple_SMS = false;
    multiple_alpha = false;
    multiple_volume = false;
    multiple_radius = false;

    // Lattice parameters
    int num_points = userInput_("lattice/numpoints", 1000);
    double first_point = userInput_("lattice/startpoint", 1.e-6);
    double lattice_size = userInput_("lattice/endpoint", 50.);
    lattice = new Lattice(num_points, first_point, lattice_size);

    //TODO: lattice = new ExpLattice(300, 1.e-5, 0.05);
}

bool Atom::Run()
{
    DebugOptions.LogFirstBuild(false);
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);

    // Check for numValenceElectrons
    numValenceElectrons_ = userInput_("NumValenceElectrons", 0);
    if(numValenceElectrons_ <= 0)
    {   *errstream << "USAGE: must have NumValenceElectrons set." << std::endl;
        return false;
    }

    // Relativistic Hartree-Fock
    core = new Core(lattice, Z, Charge);

    if(userInput_.search(2, "-c", "--clean"))
        useRead = false;
    if(userInput_.search(2, "-d", "--dont-save"))
        useWrite = false;

    if(!useRead || !ReadCore())
        core->Initialise();

    // Create Basis
    bool bsplinebasis = userInput_.search("--bspline-basis");
    bool hfbasis = userInput_.search("--hf-basis");
    bool rbasis  = userInput_.search("--r-basis");
    bool custombasis = userInput_.search("--custom-basis");

    // Generate larger mbpt basis even if this run will not actually run the MBPT part
    mbptBasisString = userInput_("MBPTBasis", "");
    bool generate_mbpt_basis = (mbptBasisString != "");

    if(bsplinebasis && !hfbasis && !rbasis && !custombasis)
        CreateBSplineBasis(generate_mbpt_basis);
    else if(!bsplinebasis && hfbasis && !rbasis && !custombasis)
        CreateHartreeFockBasis(generate_mbpt_basis);
    else if(!bsplinebasis && !hfbasis && rbasis && !custombasis)
        CreateRBasis(generate_mbpt_basis);
    else if(!bsplinebasis && !hfbasis && !rbasis && custombasis)
        CreateCustomBasis(generate_mbpt_basis);
    else
    {   *errstream << "USAGE: must have one and only one basis type (e.g. --bspline-basis)" << std::endl;
        return false;
    }

    if(useRead && useWrite)
        ReadOrWriteBasis();
    else if(useWrite)
        Write();
    else if(useRead)
        Read();

    // TODO: PrintBasis only option

    // Go do single-electron or many-electron work.
    if(numValenceElectrons_ == 1)
        RunSingleElectron();
    else
        RunMultipleElectron();

    return true;
}

void Atom::RunSingleElectron()
{
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
    NumSolutions(6), check_size_only(false), generate_mbpt_integrals(false),
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

bool Atom::ParseBasisSize(const char* basis_def, std::vector<unsigned int>& num_states)
{
    unsigned int p = 0;
    unsigned int L, max_L = 0;

    // Get Maximum L
    while(basis_def[p])
    {
        if(!isdigit(basis_def[p]))
        {
            for(L = 0; L < 10; L++)
            {   if(Constant::SpectroscopicNotation[L] == tolower(basis_def[p]))
                    break;
            }
            if(L >= 10)
                return false;
            max_L = mmax(L, max_L);
        }

        p++;
    }
    num_states.resize(max_L+1, 0);

    // Set num_states[L] to highest pqn
    p = 0;
    while(basis_def[p])
    {
        int pqn = atoi(basis_def + p);
        while(basis_def[p] && isdigit(basis_def[p]))
            p++;
        while(basis_def[p] && !isdigit(basis_def[p]))
        {
            // Get L
            for(L = 0; L < 10; L++)
            {   if(Constant::SpectroscopicNotation[L] == basis_def[p])
                    break;
            }

            num_states[L] = pqn;
            p++;
        }
    }

    // Remove core states
    std::vector<unsigned int> max_pqn_in_core(max_L+1);
    for(L = 0; L <= max_L; L++)
        max_pqn_in_core[L] = L;

    core->ToggleClosedShellCore();
    ConstStateIterator it = core->GetConstStateIterator();
    it.First();
    while(!it.AtEnd())
    {   L = it.GetStateInfo().L();
        max_pqn_in_core[L] = mmax(max_pqn_in_core[L], it.GetStateInfo().PQN());
        it.Next();
    }

    for(L = 0; L <= max_L; L++)
    {   if(max_pqn_in_core[L] > num_states[L])
            return false;

        num_states[L] -= max_pqn_in_core[L];
    }

    return true;
}

void Atom::CreateBasis(bool UseMBPT)
{
    std::vector<unsigned int> CI_basis_states;
    std::vector<unsigned int> MBPT_basis_states;

    ciBasisString = userInput_("ValenceBasis", "");
    if(ciBasisString != "")
        if(!ParseBasisSize(ciBasisString.c_str(), CI_basis_states))
        {   *errstream << "USAGE: ValenceBasis = " << ciBasisString << " has problems" << std::endl;
            exit(1);
        }

    if(UseMBPT)
    {   if(!ParseBasisSize(mbptBasisString.c_str(), MBPT_basis_states))
        {   *errstream << "USAGE: MBPTBasis = " << mbptBasisString << " has problems" << std::endl;
            exit(1);
        }

        excited_mbpt->CreateExcitedStates(MBPT_basis_states);
        *outstream << "MBPT    " << mbptBasisString << std::endl;
    }
    *outstream << "Valence " << ciBasisString << std::endl;

    excited->SetIdentifier(identifier);
    excited->CreateExcitedStates(CI_basis_states);
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
    #ifdef _MPI
        // Wait for root node to finish writing, then update integrals.
        MPI::COMM_WORLD.Barrier();
    #endif
}

bool Atom::ReadCore()
{
    // Read core electron states
    std::string filename = identifier + ".core.atom";
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        return false;

    // Check that the stored ion is the same as this one!
    double stored_Z, stored_Charge;
    fread(&stored_Z, sizeof(double), 1, fp);
    fread(&stored_Charge, sizeof(double), 1, fp);
    if((stored_Z != Z) || (stored_Charge != Charge))
    {   fclose(fp);
        *errstream << "\nAtom::ReadCore: " << filename
                   << "\n    Incorrect stored state:"
                   << "\n    Z = " << stored_Z << ", Charge = " << stored_Charge << std::endl;
        exit(1);
    }

    core->Read(fp);
    fclose(fp);
    return true;
}

bool Atom::Read()
{
    bool success = ReadCore();
    if(success && excited)
    {   // Read Excited states
        std::string filename = identifier + ".excited.atom";
        FILE* fp = fopen(filename.c_str(), "rb");
        if(!fp)
        {   return false;
        }
        excited->Read(fp);

        if(excited_mbpt)
        {   fp = freopen(filename.c_str(), "rb", fp);
            excited_mbpt->Read(fp);
        }

        fclose(fp);
    }

    return success;
}

void Atom::ReadOrWriteBasis()
{
    std::string filename = identifier + ".core.atom";
    FILE* fp = fopen(filename.c_str(), "rb");
    if(fp)
    {   fclose(fp);
        Read();
    }
    else
    {   // Wait until all nodes find there is nothing to read
        #ifdef _MPI
            MPI::COMM_WORLD.Barrier();
        #endif
        Write();
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
