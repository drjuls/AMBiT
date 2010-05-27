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
#include "Universal/Interpolator.h"
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
    multiple_length = 0;
    current_run_index = 0;
    original_id = identifier;
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

    InitialiseRunIndex();

    if(userInput_.search(2, "-c", "--clean"))
        useRead = false;

    if(userInput_.search(2, "-d", "--dont-save"))
    {   // Cannot use "-d" with multiple runs
        if(NumberRunsSelected() > 1)
        {   *errstream << "USAGE:  \"-d\" ignored: cannot use with multiple runs." << std::endl;
            useWrite = true;
        }
        else
            useWrite = false;
    }

    if(userInput_.search("HF/--read-grasp0"))
    {   // Read lattice and core and basis orbitals
        ReadGraspMCDF("MCDF.DAT");
    }
    else
    {   // Lattice parameters
        if(userInput_.search("Lattice/--exp-lattice"))
        {   int num_points = userInput_("Lattice/NumPoints", 300);
            double first_point = userInput_("Lattice/StartPoint", 1.e-5);
            double h = userInput_("Lattice/H", 0.05);
            lattice = new ExpLattice(num_points, first_point, h);
        }
        else
        {   int num_points = userInput_("Lattice/NumPoints", 1000);
            double first_point = userInput_("Lattice/StartPoint", 1.e-6);
            double lattice_size = userInput_("Lattice/EndPoint", 50.);
            lattice = new Lattice(num_points, first_point, lattice_size);
        }

        // Relativistic Hartree-Fock
        core = new Core(lattice, Z, Charge);

        InitialiseParameters();
        RunIndexBegin(true);

        if(!useRead || !ReadCore())
        {   DebugOptions.LogFirstBuild(true);
            DebugOptions.LogHFIterations(true);

            core->Initialise(userInput_("HF/Configuration", ""));
        }
    }

    // Create Basis
    bool bsplinebasis = userInput_.search("Basis/--bspline-basis");
    bool hfbasis = userInput_.search("Basis/--hf-basis");
    bool rbasis  = userInput_.search("Basis/--r-basis");
    bool custombasis = userInput_.search("Basis/--custom-basis");

    // Generate larger mbpt basis even if this run will not actually run the MBPT part
    mbptBasisString = userInput_("MBPT/Basis", "");
    bool generate_mbpt_basis = (mbptBasisString != "");

    core->ToggleOpenShellCore();

    if(bsplinebasis && !hfbasis && !rbasis && !custombasis)
        CreateBSplineBasis(generate_mbpt_basis);
    else if(!bsplinebasis && hfbasis && !rbasis && !custombasis)
        CreateHartreeFockBasis(generate_mbpt_basis);
    else if(!bsplinebasis && !hfbasis && rbasis && !custombasis)
        CreateRBasis(generate_mbpt_basis);
    else if(!bsplinebasis && !hfbasis && !rbasis && custombasis)
        CreateCustomBasis(generate_mbpt_basis);
    else
    {   *outstream << "USAGE: must have one and only one basis type (e.g. --bspline-basis)" << std::endl;
        return false;
    }

    if(useRead && useWrite)
        ReadOrWriteBasis();
    else if(useWrite)
        Write();
    else if(useRead)
        Read();

    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(false);

    // Loop through all multiple runs saving the basis
    if(NumberRunsSelected() > 1)
    {
        RunIndexNext(false);
        while(!RunIndexAtEnd())
        {
            core->ToggleOpenShellCore();

            // Try read, if it fails then update and write
            if(!useRead || !Read())
            {
                if(!useRead || !ReadCore())
                    core->Update();
                excited->Update();
                if(excited_mbpt)
                    excited_mbpt->Update();
                Write();
            }

            RunIndexNext(false);
        }
        RunIndexBegin(false);

        // From now on, read core and basis
        useRead = true;
    }

    // Print basis only option
    if(userInput_.search(2, "--print-basis", "-p"))
    {   // Check follower for option
        std::string print_option = userInput_.next("");
        if(print_option == "Cowan")
            GenerateCowanInputFile();
        else
            WriteGraspMCDF();

        return true;
    }

    // Go do single-electron or many-electron work.
    if(numValenceElectrons_ == 1)
        RunSingleElectron();
    else
        RunMultipleElectron();

    return true;
}

void Atom::RunSingleElectron()
{
    core->ToggleClosedShellCore();

    mbpt = NULL;
    if(excited_mbpt)
        mbpt = new CoreMBPTCalculator(lattice, core, excited_mbpt);

    bool sigma_potential = userInput_.search("Basis/Valence/--sigma-potential") && mbpt;

    while(!RunIndexAtEnd())
    {
        SetRunParameters(true);
        SetRunCore();

        if(mbpt)
            mbpt->UpdateIntegrals(excited);

        StateIterator it = excited->GetStateIterator();
        it.First();
        while(!it.AtEnd())
        {
            StateInfo si = it.GetStateInfo();
            DiscreteState* ds = it.GetState();
            *outstream << "\n" << ds->Name() << "\n";
            *outstream << std::setprecision(12);
            *outstream << "HF Energy: " << ds->Energy() * Constant::HartreeEnergy_cm << std::endl;

            if(sigma_potential)
            {
                // Create sigma, if it doesn't exist, and find matrix element
                double dE = 0.0;
                if(!excited->RetrieveSecondOrderSigma(si))
                    dE = excited->CreateSecondOrderSigma(si, *mbpt);
                else
                    dE = excited->GetSigmaMatrixElement(si);
                *outstream << "HF + MBPT: " << (ds->Energy() + dE) * Constant::HartreeEnergy_cm << std::endl;        
                
                // Iterate to create Brueckner orbital
                DiscreteState brueckner = excited->GetStateWithSigma(si);
                *outstream << "Brueckner: " << brueckner.Energy() * Constant::HartreeEnergy_cm << std::endl;

                // Set energy to known value, if requested
                double E_experiment = userInput_(("Basis/Valence/" + si.Name()).c_str(), 0.0);
                if(E_experiment)
                {   excited->SetEnergyViaSigma(si, E_experiment/Constant::HartreeEnergy_cm);
                    brueckner = excited->GetStateWithSigma(si);
                }

                // Log ratio
                *logstream << "\n" << brueckner.Name() << std::endl;
                for(unsigned int i = 0; i < 100; i+=10)
                {
                    *logstream << lattice->R(i) << "\t";
                    *logstream << brueckner.f[i]/ds->f[i] << "\t" << brueckner.g[i]/ds->g[i] << std::endl;
                }
            }
            else if(mbpt)
            {   double dE = mbpt->GetOneElectronDiagrams(si, si);
                *outstream << "HF + MBPT: " << (ds->Energy() + dE) * Constant::HartreeEnergy_cm << std::endl;        
            }

            it.Next();
        }

        excited->ClearSigmas();
        RunIndexNext(false);
    }
    RunIndexBegin(false);
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

    ciBasisString = userInput_("Basis/ValenceBasis", "");
    if(ciBasisString != "")
        if(!ParseBasisSize(ciBasisString.c_str(), CI_basis_states))
        {   *errstream << "USAGE: Basis/ValenceBasis = " << ciBasisString << " has problems" << std::endl;
            exit(1);
        }

    if(UseMBPT)
    {   if(!ParseBasisSize(mbptBasisString.c_str(), MBPT_basis_states))
        {   *errstream << "USAGE: MBPT/Basis = " << mbptBasisString << " has problems" << std::endl;
            exit(1);
        }

        core->ToggleOpenShellCore();
        excited_mbpt->CreateExcitedStates(MBPT_basis_states);
        *outstream << "MBPT    " << mbptBasisString << std::endl;
    }
    *outstream << "Valence " << ciBasisString << std::endl;

    excited->SetIdentifier(identifier);
    core->ToggleOpenShellCore();
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

    unsigned int n = userInput_("Basis/BSpline/N", 40);
    unsigned int k = userInput_("Basis/BSpline/K", 7);
    double Rmax = userInput_("Basis/BSpline/Rmax", 50.);
    basis->SetParameters(n, k, Rmax);
    CreateBasis(UseMBPT);
}

void Atom::CreateCustomBasis(bool UseMBPT)
{
    std::string filename = userInput_("Basis/CustomFile", "CustomBasis.txt");
    CustomBasis* basis;

    if(UseMBPT)
    {   excited_mbpt = new CustomBasis(lattice, core);
        excited = new SubsetBasis(lattice, excited_mbpt);
        basis = dynamic_cast<CustomBasis*>(excited_mbpt);
    }
    else
    {   excited = new CustomBasis(lattice, core);
        basis = dynamic_cast<CustomBasis*>(excited);
    }

    basis->SetFile(filename);
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

void Atom::CreateReadBasis(bool UseMBPT)
{}

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
    FILE* fp = fopen("mitch.txt", "wt");

    const DiscreteState* ds = core->GetState(StateInfo(1, -1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(2, -1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(2, 1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(StateInfo(2, -2));
    PrintWavefunctionCowan(fp, ds);

    for(int n = 3; n <= 5; n++)
    {
        ds = excited->GetState(StateInfo(n, -1));
        PrintWavefunctionCowan(fp, ds);
//        ds = excited->GetState(StateInfo(n, 1));
//        PrintWavefunctionCowan(fp, ds);
//        ds = excited->GetState(StateInfo(n, -2));
//        PrintWavefunctionCowan(fp, ds);
        ds = excited->GetState(StateInfo(n, 2));
        PrintWavefunctionCowan(fp, ds);
        ds = excited->GetState(StateInfo(n, -3));
        PrintWavefunctionCowan(fp, ds);
//        if(n > 3)
//        {
//            ds = excited->GetState(StateInfo(n, 3));
//            PrintWavefunctionCowan(fp, ds);
//            ds = excited->GetState(StateInfo(n, -4));
//            PrintWavefunctionCowan(fp, ds);
//        }
    }

    return;
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
/*
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
*/
    // Upper component
    fprintf(fp, "     %2d%2d     %5d\n", ds->RequiredPQN(), ds->Kappa(), ds->Size());

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
/*    while(i < lattice->Size())
    {   if(count == 5)
        {   fprintf(fp, "\n");
            count = 1;
        }
        else
            count++;
        fprintf(fp, "%14.7E", 0.0);
        i++;
    }
/*    
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
*/
    //fprintf(fp, "\n     %4s     %5d%12.4E%12.4E\n", ds->Name().c_str(), lattice->Size(), f[0], f[1]);
    fprintf(fp, "\n     %2d%2d     %5d\n", ds->RequiredPQN(), ds->Kappa(), ds->Size());

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
/*    while(i < lattice->Size())
    {   if(count == 5)
        {   fprintf(fp, "\n");
            count = 1;
        }
        else
            count++;
        fprintf(fp, "%14.7E", 0.0);
        i++;
    }
*/    
    fprintf(fp, "\n");
}

bool Atom::ReadGraspMCDF(const std::string& filename)
{
    FILE* fp = fopen(filename.c_str(), "rb");
    int record_size;    // Need to read size of record (bytes) at start and end of each record
    unsigned int i;

    // Read Record 1 (see GRASP0 documentation from Patrick Norrington's website).
    char IHED[80], RECORD[20];
    fread(&record_size, sizeof(int), 1, fp);
    fread(IHED, sizeof(char), 80, fp);
    fread(RECORD, sizeof(char), 20, fp);
    fread(&record_size, sizeof(int), 1, fp);
    i = 80;
    while(i>0 && isspace(IHED[--i]))
        IHED[i] = 0;

    // Record 2
    int NCMIN;
    int NW;     // total number of states
    int NCF;
    int N;      // lattice size
    fread(&record_size, sizeof(int), 1, fp);
    fread(&NCMIN, sizeof(int), 1, fp);
    fread(&NW, sizeof(int), 1, fp);
    fread(&NCF, sizeof(int), 1, fp);
    fread(&N, sizeof(int), 1, fp);
    fread(&record_size, sizeof(int), 1, fp);

    // Record 3
    double RNT; // lattice->R(0)
    double H;   // lattice->H()
    double C;   // 1/Constant::Alpha
    fread(&record_size, sizeof(int), 1, fp);
    fread(&Z, sizeof(double), 1, fp);
    fread(&RNT, sizeof(double), 1, fp);
    fread(&H, sizeof(double), 1, fp);
    fread(&C, sizeof(double), 1, fp);
    fread(&record_size, sizeof(int), 1, fp);

    // Create lattice and core
    lattice = new ExpLattice(N, RNT, H);
    Constant::Alpha = 1./C;
    Constant::AlphaSquared = Constant::Alpha*Constant::Alpha;
    *outstream << "Speed of light (1/alpha) = " << C << std::endl;
    core = new Core(lattice, Z, Charge);

    InitialiseParameters();
    RunIndexBegin();

    // Create Basis
    if(userInput_.search("Basis/--hf-basis"))
        excited = new HartreeFockBasis(lattice, core);
    else if(userInput_.search("Basis/--r-basis"))
        excited = new RStates(lattice, core);
    else if(userInput_.search("Basis/--custom-basis"))
        excited = new CustomBasis(lattice, core);
    else
    {   *outstream << "USAGE: --read-grasp0 requires basis specified as\n"
                   << "       --hf-basis, --r-basis, or --custom-basis." << std::endl;
        exit(1);
    }

    // Record 4 - 6, repeated for each orbital
    unsigned int orbital_count = 0;

    // Need to know where the orbital belongs: core or valence or excited
    std::string configuration = userInput_("HF/Configuration", "");
    if(!configuration.size())
    {   *outstream << "USAGE: --read-grasp0 requires \"HFConfiguration\" set in input file." << std::endl;
        exit(1);
    }
    std::string closed_shell_string;
    std::string open_shell_string;

    unsigned int colon_pos = configuration.find(':');
    if(colon_pos == std::string::npos)
        closed_shell_string = configuration;
    else
    {   closed_shell_string = configuration.substr(0, colon_pos);
        open_shell_string = configuration.substr(colon_pos+1, configuration.size()-colon_pos-1);
    }

    Configuration closed_shell_config(closed_shell_string);
    Configuration open_shell_config(open_shell_string);

    if(closed_shell_config.NumParticles() + open_shell_config.NumParticles() != Z - Charge)
    {   *errstream << "Core::BuildFirstApproximation: Incorrect electron count in configuration." << std::endl;
        exit(1);
    }
    
    while(orbital_count < NW)
    {
        // Record 4
        char NH[2];
        int NP;     // ds->RequiredPQN()
        int NAK;    // ds->Kappa()
        double E;   // fabs(ds->Energy())
        fread(&record_size, sizeof(int), 1, fp);
        fread(&NH, sizeof(char), 2, fp);
        fread(&NP, sizeof(int), 1, fp);
        fread(&NAK, sizeof(int), 1, fp);
        fread(&E, sizeof(double), 1, fp);
        fread(&record_size, sizeof(int), 1, fp);

        DiscreteState* ds = new DiscreteState(NP, NAK);
        ds->SetEnergy(-E);

        // Record 5
        double PZ;
        double QZ;
        record_size = 2*sizeof(double);
        fread(&record_size, sizeof(int), 1, fp);
        fread(&PZ, sizeof(double), 1, fp);
        fread(&QZ, sizeof(double), 1, fp);
        fread(&record_size, sizeof(int), 1, fp);        

        // Record 6
        double P[N], Q[N];  // upper and lower components
        fread(&record_size, sizeof(int), 1, fp);
        fread(P, sizeof(double), N, fp);
        fread(Q, sizeof(double), N, fp);
        fread(&record_size, sizeof(int), 1, fp);
        
        i = N - 1;
        while((i > 0) && (P[i] == 0.0) && (Q[i] == 0.0))
            i--;
        ds->ReSize(i+1);
        for(i = 0; i < ds->Size(); i++)
        {   ds->f[i] = P[i];
            ds->g[i] = Q[i]/Constant::Alpha;
        }

        // Get derivatives
        Interpolator interp(lattice);
        unsigned int order = 6;
        interp.GetDerivative(ds->f, ds->df, order);
        interp.GetDerivative(ds->g, ds->dg, order);

        // Split electrons between subshells
        double occupancy = 2.*fabs(ds->Kappa()) / (4.*ds->L()+2.);

        // Check non-rel configurations with a non-rel info
        NonRelInfo non_rel_info(ds->RequiredPQN(), ds->L());

        if(closed_shell_config.GetOccupancy(non_rel_info))
        {   ds->SetOccupancy(occupancy * closed_shell_config.GetOccupancy(non_rel_info));
            core->AddState(ds);
        }
        else
        {   if(open_shell_config.GetOccupancy(non_rel_info))
            {   ds->SetOccupancy(occupancy * open_shell_config.GetOccupancy(non_rel_info));
                core->AddState(ds);
                core->SetOpenShellState(ds, ds->Occupancy());
            }

            excited->AddState(ds);
        }
        *outstream << "Read state " << ds->Name() << " " << ds->Energy() << std::endl;

        orbital_count++;
    }

    fclose(fp);
    return true;
}

void Atom::WriteGraspMCDF() const
{
    FILE* fp = fopen("MCDF.DMP", "wb");
    int record_size;    // Need to write size of record (bytes) at start and end of each record
    
    // Write Record 1 (see GRASP0 documentation from Patrick Norrington's website).
    record_size = 100;
    char IHED[80], RECORD[20];
    memset(IHED, ' ', 80);
    memset(RECORD, ' ', 20);
    memcpy(IHED, identifier.c_str(), identifier.size());
    fwrite(&record_size, sizeof(int), 1, fp);
    fwrite(IHED, sizeof(char), 80, fp);
    fwrite(RECORD, sizeof(char), 20, fp);
    fwrite(&record_size, sizeof(int), 1, fp);
    
    // Record 2
    core->ToggleClosedShellCore();  // Avoid double counting valence states
    int NCMIN = 0;
    int NW = core->NumStates() + excited->NumStates();
    int NCF = 0;
    int N = 0;      // Largest orbital size
    ConstStateIterator it = core->GetConstStateIterator();
    it.First();
    while(!it.AtEnd())
    {   N = mmax(N, it.GetState()->Size());
        it.Next();
    }
    it = excited->GetConstStateIterator();
    it.First();
    while(!it.AtEnd())
    {   N = mmax(N, it.GetState()->Size());
        it.Next();
    }
    record_size = 4 * sizeof(int);
    fwrite(&record_size, sizeof(int), 1, fp);
    fwrite(&NCMIN, sizeof(int), 1, fp);
    fwrite(&NW, sizeof(int), 1, fp);
    fwrite(&NCF, sizeof(int), 1, fp);
    fwrite(&N, sizeof(int), 1, fp);
    fwrite(&record_size, sizeof(int), 1, fp);
    
    // Record 3
    double RNT = lattice->R(0);
    double H = lattice->H();
    double C = 1/Constant::Alpha;
    record_size = 4*sizeof(double);
    fwrite(&record_size, sizeof(int), 1, fp);
    fwrite(&Z, sizeof(double), 1, fp);
    fwrite(&RNT, sizeof(double), 1, fp);
    fwrite(&H, sizeof(double), 1, fp);
    fwrite(&C, sizeof(double), 1, fp);
    fwrite(&record_size, sizeof(int), 1, fp);

    *logstream << "Writing grasp file. Num orbitals = " << NW << std::endl;
    // Record 4 - 6, repeated for each orbital
    it = core->GetConstStateIterator();
    it.First();
    while(!it.AtEnd())
    {
        const DiscreteState* ds = it.GetState();
        *logstream << " " << ds->Name() << std::endl;
        WriteGraspMcdfOrbital(fp, ds, N);
        it.Next();
    }

    it = excited->GetConstStateIterator();
    it.First();
    while(!it.AtEnd())
    {
        const DiscreteState* ds = it.GetState();
        *logstream << " " << ds->Name() << std::endl;
        WriteGraspMcdfOrbital(fp, ds, N);
        it.Next();
    }

    fclose(fp);
}

void Atom::WriteGraspMcdfOrbital(FILE* fp, const DiscreteState* ds, unsigned int lattice_size) const
{
    unsigned int record_size;
    unsigned int N = lattice_size;

    // Record 4
    char NH[2];
    int NP = ds->RequiredPQN();
    int NAK = ds->Kappa();
    double E = fabs(ds->Energy());
    NH[0] = toupper(Constant::SpectroscopicNotation[ds->L()]);
    if(NAK > 0)
        NH[1] = '-';
    else
        NH[1] = ' ';
    record_size = 2 + 2*sizeof(int) + sizeof(double);
    fwrite(&record_size, sizeof(int), 1, fp);
    fwrite(&NH, sizeof(char), 2, fp);
    fwrite(&NP, sizeof(int), 1, fp);
    fwrite(&NAK, sizeof(int), 1, fp);
    fwrite(&E, sizeof(double), 1, fp);
    fwrite(&record_size, sizeof(int), 1, fp);

    // Copy upper and lower components to buffers P and Q
    double P[N], Q[N];
    bool switch_sign = (ds->f[0] < 0.0);
    unsigned int i;
    for(i = 0; i < ds->Size(); i++)
    {   if(switch_sign)
        {   P[i] = -ds->f[i];
            Q[i] = -ds->g[i]*Constant::Alpha;
        }
        else
        {   P[i] = ds->f[i];
            Q[i] = ds->g[i]*Constant::Alpha;
        }
    }
    while(i < N)
    {   P[i] = 0.0;
        Q[i] = 0.0;
        i++;
    }

    // Record 5
    double PZ = P[0]/lattice->R(0);
    double QZ = Q[0]/lattice->R(0);
    record_size = 2*sizeof(double);
    fwrite(&record_size, sizeof(int), 1, fp);
    fwrite(&PZ, sizeof(double), 1, fp);
    fwrite(&QZ, sizeof(double), 1, fp);
    fwrite(&record_size, sizeof(int), 1, fp);        
   
    // Record 6
    record_size = 2*N*sizeof(double);
    fwrite(&record_size, sizeof(int), 1, fp);
    fwrite(P, sizeof(double), N, fp);
    fwrite(Q, sizeof(double), N, fp);
    fwrite(&record_size, sizeof(int), 1, fp);        
}
