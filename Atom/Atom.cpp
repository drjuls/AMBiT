#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Atom.h"
#include "Basis/BasisGenerator.h"

#include "MBPT/CoreMBPTCalculator.h"
#include "Universal/ScalapackMatrix.h"
#include "Universal/ExpLattice.h"
#include "Universal/Eigensolver.h"
#include "Universal/Interpolator.h"
#include "RateCalculator.h"
#include <fstream>
#include <sstream>

Atom::Atom(const MultirunOptions userInput, unsigned int atomic_number, const std::string& atom_identifier):
    user_input(userInput), Z(atomic_number), identifier(atom_identifier)
{
    lattice = pLattice();
    integrals = nullptr;
    orbitals = nullptr;

//    integralsMBPT = NULL;
//    mbpt = NULL;
//    valence_mbpt = NULL;
//    sigma3 = NULL;

    // Operational parameters
    useRead = true;
    useWrite = true;

    // CI + MBPT parameters
    NumSolutions = 6;
    MaxEnergy = 0.0;
    mbptBasisString = "";
    ciBasisString = "";
    check_size_only = false;
    save_eigenstates = true;
    generate_mbpt_integrals = false;
    includeSigma1 = false;
    includeSigma2 = false;
    includeSigma3 = false;
}

Atom::~Atom(void)
{
    if(integrals)
        delete integrals;
}

void Atom::CreateBasis()
{
    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);

    int N = user_input("HF/N", -1);  // Number of electrons
    if(N == -1 || N > Z)
    {   *errstream << "Must specify N (number of electrons) with Z >= N" << std::endl;
        exit(1);
    }

    // Check for numValenceElectrons
    num_valence_electrons = user_input("NumValenceElectrons", 0);
    if(num_valence_electrons <= 0)
    {   *errstream << "USAGE: must have NumValenceElectrons set." << std::endl;
        exit(1);
    }

    if(user_input.search(2, "-c", "--clean"))
        useRead = false;

    if(user_input.search(2, "-d", "--dont-save"))
    {   // Cannot use "-d" with multiple runs
        if(NumberRunsSelected() > 1)
        {   *errstream << "USAGE:  \"-d\" ignored: cannot use with multiple runs." << std::endl;
            useWrite = true;
        }
        else
            useWrite = false;
    }

    if(user_input.search("HF/--read-grasp0"))
    {   // Read lattice and core and basis orbitals
        ReadGraspMCDF("MCDF.DAT");
    }
    else
    {   // Lattice parameters
        if(user_input.search("Lattice/--exp-lattice"))
        {   int num_points = user_input("Lattice/NumPoints", 300);
            double first_point = user_input("Lattice/StartPoint", 1.e-5);
            double h = user_input("Lattice/H", 0.05);
            lattice = pLattice(new ExpLattice(num_points, first_point, h));
        }
        else
        {   int num_points = user_input("Lattice/NumPoints", 1000);
            double first_point = user_input("Lattice/StartPoint", 1.e-6);
            double lattice_size = user_input("Lattice/EndPoint", 50.);
            lattice = pLattice(new Lattice(num_points, first_point, lattice_size));
        }
    }

    // Relativistic Hartree-Fock
    BasisGenerator basis_generator(lattice, user_input);
    hf_core = basis_generator.GenerateHFCore();

    // Create Basis
    orbitals = basis_generator.GenerateBasis();

    if(useRead && useWrite)
        ReadOrWriteBasis();
    else if(useWrite)
        Write();
    else if(useRead)
        Read();
}

/*
void Atom::RunSingleElectron()
{
    core->ToggleClosedShellCore();

    mbpt = NULL;
    if(excited_mbpt)
        mbpt = new CoreMBPTCalculator(lattice, core, excited_mbpt);
    
    MathConstant* constants = MathConstant::Instance();

    bool sigma_potential = userInput_.search("Basis/Valence/--sigma-potential") && mbpt;

    while(!RunIndexAtEnd())
    {
        SetRunParameters(true);
        SetRunCore();

        if(mbpt)
            mbpt->UpdateIntegrals(excited);

        StateIterator it = excited->GetStateIterator();
        it.First();
        int ZIterations = (int) userInput_("ZIterations", 1);
        bool UseRecursive = userInput_.search("--recursive-build");
        if(ZIterations == 1 && !UseRecursive)
        {
            while(!it.AtEnd())
            {
                OrbitalInfo si = it.GetOrbitalInfo();
                pOrbital ds = it.GetState();
                *outstream << "\n" << ds->Name() << "\n";
                *outstream << std::setprecision(12);
                *outstream << "HF Energy: " << ds->Energy() * constants->HartreeEnergyInInvCm() << std::endl;
    
                if(sigma_potential)
                {
                    // Create sigma, if it doesn't exist, and find matrix element
                    double dE = 0.0;
                    if(!excited->RetrieveSecondOrderSigma(si))
                        dE = excited->CreateSecondOrderSigma(si, *mbpt);
                    else
                        dE = excited->GetSigmaMatrixElement(si);
                    *outstream << "HF + MBPT: " << (ds->Energy() + dE) * constants->HartreeEnergyInInvCm() << std::endl;
                    
                    // Iterate to create Brueckner orbital
                    Orbital brueckner = excited->GetStateWithSigma(si);
                    *outstream << "Brueckner: " << brueckner.Energy() * constants->HartreeEnergyInInvCm() << std::endl;
    
                    // Set energy to known value, if requested
                    double E_experiment = userInput_(("Basis/Valence/" + si.Name()).c_str(), 0.0);
                    if(E_experiment)
                    {   excited->SetEnergyViaSigma(si, E_experiment/constants->HartreeEnergyInInvCm());
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
                    *outstream << "HF + MBPT: " << (ds->Energy() + dE) * constants->HartreeEnergyInInvCm() << std::endl;
                }
    
                it.Next();
            }
        }
        else if(ZIterations > 1)
        {
            if(userInput_("Z", Z) == Z)
            {
                *outstream << std::endl;
                *outstream << identifier << " ";

                it = core->GetStateIterator();
                it.First();
                while(!it.AtEnd())
                {
                    *outstream << it.GetState()->Name() << " ";
                    it.Next();
                }
                it = excited->GetStateIterator();
                it.First();
                while(!it.AtEnd())
                {
                    *outstream << it.GetState()->Name() << " ";
                    it.Next();
                }
                *outstream << std::endl;
            }
            
            *outstream << Z << " ";
            
            it = core->GetStateIterator();
            it.First();

            while(!it.AtEnd())
            {
                OrbitalInfo si = it.GetOrbitalInfo();
                pOrbital ds = it.GetState();
                *outstream << std::setprecision(12);

                if(mbpt)
                {   double dE = mbpt->GetOneElectronDiagrams(si, si);
                    *outstream << (ds->Energy() + dE) * constants->HartreeEnergyInInvCm() << " ";
                }
                else
                {
                    *outstream << ds->Energy() * constants->HartreeEnergyInInvCm() << " ";
                }
    
                it.Next();
            }

            it = excited->GetStateIterator();
            it.First();

            while(!it.AtEnd())
            {
                OrbitalInfo si = it.GetOrbitalInfo();
                pOrbital ds = it.GetState();
                *outstream << std::setprecision(12);

                if(mbpt)
                {   double dE = mbpt->GetOneElectronDiagrams(si, si);
                    *outstream << (ds->Energy() + dE) * constants->HartreeEnergyInInvCm() << " ";
                }
                else
                {
                    *outstream << ds->Energy() * constants->HartreeEnergyInInvCm() << " ";
                }
    
                it.Next();
            }

            *outstream << std::endl;
        }
        else if(UseRecursive)
        {
            if(true)
            {
                *outstream << identifier << " ";

                it = core->GetStateIterator();
                it.First();
                while(!it.AtEnd())
                {
                    *outstream << it.GetState()->Name() << " ";
                    it.Next();
                }
                it = excited->GetStateIterator();
                it.First();
                while(!it.AtEnd())
                {
                    *outstream << it.GetState()->Name() << " ";
                    it.Next();
                }
                *outstream << std::endl;
            }
            
            *outstream << userInput_("N", 1) << " ";
            
            it = core->GetStateIterator();
            it.First();

            while(!it.AtEnd())
            {
                OrbitalInfo si = it.GetOrbitalInfo();
                pOrbital ds = it.GetState();
                *outstream << std::setprecision(12);

                if(mbpt)
                {   double dE = mbpt->GetOneElectronDiagrams(si, si);
                    *outstream << (ds->Energy() + dE) * constants->HartreeEnergyInInvCm() << " ";
                }
                else
                {
                    *outstream << ds->Energy() * constants->HartreeEnergyInInvCm() << " ";
                }
    
                it.Next();
            }

            it = excited->GetStateIterator();
            it.First();

            while(!it.AtEnd())
            {
                OrbitalInfo si = it.GetOrbitalInfo();
                pOrbital ds = it.GetState();
                *outstream << std::setprecision(12);

                if(mbpt)
                {   double dE = mbpt->GetOneElectronDiagrams(si, si);
                    *outstream << (ds->Energy() + dE) * constants->HartreeEnergyInInvCm() << " ";
                }
                else
                {
                    *outstream << ds->Energy() * constants->HartreeEnergyInInvCm() << " ";
                }
    
                it.Next();
            }

            *outstream << std::endl;
        }

        excited->ClearSigmas();
        RunIndexNext(false);
    }
    RunIndexBegin(false);
}

StateIterator Atom::GetIteratorToNextOrbitalToFill()
{
    core->ToggleClosedShellCore();
    StateIterator si = core->GetStateIterator();
    si.First();

    std::string result = "";
    double lowestEnergy = 0;
    StateIterator lowestsi = si;

    while(!si.AtEnd())
    {
        OrbitalInfo oi = si.GetOrbitalInfo();
        pOrbital orb = si.GetState();
        if((orb->Energy() < lowestEnergy) && ((int) orb->Occupancy() < (int) oi.MaxNumElectrons()))
        {
            lowestEnergy = orb->Energy();
            lowestsi = si;
        }
        si.Next();
    }
    si = excited->GetStateIterator();
    si.First();
    while(!si.AtEnd())
    {
        OrbitalInfo oi = si.GetOrbitalInfo();
        pOrbital orb = si.GetState();
        if((orb->Energy() < lowestEnergy))
        {
            lowestEnergy = orb->Energy();
            lowestsi = si;
        }
        si.Next();
    }

    return lowestsi;
}

std::string Atom::GetNextConfigString()
{
    std::stringstream configstring (std::stringstream::in | std::stringstream::out);
    std::stringstream resultstring (std::stringstream::in | std::stringstream::out);
    configstring << userInput_("HF/Configuration", "");
    int mPQN, aPQN;
    char mL, aL;
    int mOccupancy;
    bool foundExisting = false;

    std::stringstream addthis (std::stringstream::in | std::stringstream::out);
    addthis << GetIteratorToNextOrbitalToFill().GetState()->Name();
    addthis >> aPQN >> aL;
    
    configstring.seekg(0, std::ios::end);
    int length = configstring.tellg();
    configstring.seekg(0, std::ios::beg);
    
    while(configstring.tellg() < length && (int) configstring.tellg() != (int) -1)
    {
        mPQN = 0;
        mOccupancy = 0;
        configstring >> mPQN >> mL >> mOccupancy;
        if(mPQN == aPQN && mL == aL)
        {
            resultstring << mPQN << mL << mOccupancy + 1 << " ";
            foundExisting = true;
        }
        else
        {
            if(mOccupancy > 0)
            {
                resultstring << mPQN << mL << mOccupancy << " ";
            }
        }
    }
    if(foundExisting == false)
    {
        resultstring << aPQN << aL << 1 << " ";
    }

    return resultstring.str();
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
            L = MathConstant::Instance()->GetL(tolower(basis_def[p]));
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
            L = MathConstant::Instance()->GetL(tolower(basis_def[p]));
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
    {   L = it.GetOrbitalInfo().L();
        max_pqn_in_core[L] = mmax(max_pqn_in_core[L], it.GetOrbitalInfo().PQN());
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
    if(userInput_("Z", Z) == Z && !userInput_.search("--recursive-build"))
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

    BSplineBasis::SplineType spline_type = BSplineBasis::Reno;
    std::string spline = userInput_("Basis/BSpline/SplineType", "Reno");
    if(spline.compare("Reno") == 0 || spline.compare("DKB") == 0)
        spline_type = BSplineBasis::Reno;
    else if(spline.compare("Vanderbilt") == 0)
        spline_type = BSplineBasis::Vanderbilt;
    else if(spline.compare("NotreDame") == 0 || spline.compare("Johnson") == 0)
        spline_type = BSplineBasis::NotreDame;

    basis->SetParameters(n, k, Rmax, spline_type);

    bool reorth = userInput_.search("Basis/BSpline/--reorthogonalise");
    basis->OrthogonaliseAgain(reorth);

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
        excited->SetIdentifier(identifier);

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
*/
void Atom::GenerateCowanInputFile()
{
    FILE* fp = fopen("mitch.txt", "wt");

    pOrbitalConst ds = core->GetState(OrbitalInfo(1, -1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(OrbitalInfo(2, -1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(OrbitalInfo(2, 1));
    PrintWavefunctionCowan(fp, ds);
    ds = excited->GetState(OrbitalInfo(2, -2));
    PrintWavefunctionCowan(fp, ds);

    for(int n = 3; n <= 5; n++)
    {
        ds = excited->GetState(OrbitalInfo(n, -1));
        PrintWavefunctionCowan(fp, ds);
//        ds = excited->GetState(OrbitalInfo(n, 1));
//        PrintWavefunctionCowan(fp, ds);
//        ds = excited->GetState(OrbitalInfo(n, -2));
//        PrintWavefunctionCowan(fp, ds);
        ds = excited->GetState(OrbitalInfo(n, 2));
        PrintWavefunctionCowan(fp, ds);
        ds = excited->GetState(OrbitalInfo(n, -3));
        PrintWavefunctionCowan(fp, ds);
//        if(n > 3)
//        {
//            ds = excited->GetState(OrbitalInfo(n, 3));
//            PrintWavefunctionCowan(fp, ds);
//            ds = excited->GetState(OrbitalInfo(n, -4));
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

void Atom::PrintWavefunctionCowan(FILE* fp, pOrbitalConst ds)
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
    fprintf(fp, "     %4s     %5d%12.4E%12.4E%12.4E\n", ds->Name().c_str(), lattice->size(), f[0], f[1], -ds->Energy());
*/
    // Upper component
    fprintf(fp, "     %2d%2d     %5d\n", ds->PQN(), ds->Kappa(), ds->size());

    unsigned int count = 0;
    unsigned int i;
    for(i = 0; i < ds->size(); i++)
    {   if(count == 5)
        {   fprintf(fp, "\n");
            count = 1;
        }
        else
            count++;
        fprintf(fp, "%14.7E", ds->f[i]);
    }
/*    while(i < lattice->size())
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
    //fprintf(fp, "\n     %4s     %5d%12.4E%12.4E\n", ds->Name().c_str(), lattice->size(), f[0], f[1]);
    fprintf(fp, "\n     %2d%2d     %5d\n", ds->PQN(), ds->Kappa(), ds->size());

    count = 0;
    for(i = 0; i < ds->size(); i++)
    {   if(count == 5)
        {   fprintf(fp, "\n");
            count = 1;
        }
        else
            count++;
        fprintf(fp, "%14.7E", ds->g[i]*PhysicalConstant::Instance()->GetAlpha());
    }
/*    while(i < lattice->size())
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
    double C;   // 1/alpha
    fread(&record_size, sizeof(int), 1, fp);
    fread(&Z, sizeof(double), 1, fp);
    fread(&RNT, sizeof(double), 1, fp);
    fread(&H, sizeof(double), 1, fp);
    fread(&C, sizeof(double), 1, fp);
    fread(&record_size, sizeof(int), 1, fp);

    // Create lattice and core
    lattice = pLattice(new ExpLattice(N, RNT, H));
    PhysicalConstant::Instance()->SetAlpha(1./C);
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

    size_t colon_pos = configuration.find(':');
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
        int NP;     // ds->PQN()
        int NAK;    // ds->Kappa()
        double E;   // fabs(ds->Energy())
        fread(&record_size, sizeof(int), 1, fp);
        fread(&NH, sizeof(char), 2, fp);
        fread(&NP, sizeof(int), 1, fp);
        fread(&NAK, sizeof(int), 1, fp);
        fread(&E, sizeof(double), 1, fp);
        fread(&record_size, sizeof(int), 1, fp);

        pOrbital ds(new Orbital(NAK, NP, -E));
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
        ds->size(i+1);
        for(i = 0; i < ds->size(); i++)
        {   ds->f[i] = P[i];
            ds->g[i] = Q[i]/PhysicalConstant::Instance()->GetAlpha();
        }

        // Get derivatives
        Interpolator interp(lattice);
        unsigned int order = 6;
        interp.GetDerivative(ds->f, ds->dfdr, order);
        interp.GetDerivative(ds->g, ds->dgdr, order);

        // Split electrons between subshells
        double occupancy = 2.*fabs(ds->Kappa()) / (4.*ds->L()+2.);

        // Check non-rel configurations with a non-rel info
        NonRelInfo non_rel_info(ds->PQN(), ds->L());

        if(closed_shell_config.GetOccupancy(non_rel_info))
        {   ds->SetOccupancy(occupancy * closed_shell_config.GetOccupancy(non_rel_info));
            core->AddState(ds);
        }
        else
        {   if(open_shell_config.GetOccupancy(non_rel_info))
            {   ds->SetOccupancy(occupancy * open_shell_config.GetOccupancy(non_rel_info));
                core->AddState(ds);
                core->SetOpenShellState(OrbitalInfo(ds), ds->Occupancy());
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
    {   N = mmax(N, it.GetState()->size());
        it.Next();
    }
    it = excited->GetConstStateIterator();
    it.First();
    while(!it.AtEnd())
    {   N = mmax(N, it.GetState()->size());
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
    double C = 1/PhysicalConstant::Instance()->GetAlpha();
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
        pOrbitalConst ds = it.GetState();
        *logstream << " " << ds->Name() << std::endl;
        WriteGraspMcdfOrbital(fp, ds, N);
        it.Next();
    }

    it = excited->GetConstStateIterator();
    it.First();
    while(!it.AtEnd())
    {
        pOrbitalConst ds = it.GetState();
        *logstream << " " << ds->Name() << std::endl;
        WriteGraspMcdfOrbital(fp, ds, N);
        it.Next();
    }

    fclose(fp);
}

void Atom::WriteGraspMcdfOrbital(FILE* fp, pOrbitalConst ds, unsigned int lattice_size) const
{
    unsigned int record_size;
    unsigned int N = lattice_size;

    // Record 4
    char NH[2];
    int NP = ds->PQN();
    int NAK = ds->Kappa();
    double E = fabs(ds->Energy());
    NH[0] = toupper(MathConstant::Instance()->GetSpectroscopicNotation(ds->L()));
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
    for(i = 0; i < ds->size(); i++)
    {   if(switch_sign)
        {   P[i] = -ds->f[i];
            Q[i] = -ds->g[i]*PhysicalConstant::Instance()->GetAlpha();
        }
        else
        {   P[i] = ds->f[i];
            Q[i] = ds->g[i]*PhysicalConstant::Instance()->GetAlpha();
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

void Atom::WriteEigenstatesToSolutionMap()
{
    if(GetSymmetryEigenstatesMap()->empty())
    {
        *errstream << "Could not find eigenstates. Generate or read them first." << std::endl;
        exit(1);
    }
    if(!GetSolutionMap()->empty())
    {
        GetSolutionMap()->clear();
    }
    if(!GetSymmetryEigenstatesMap()->empty() && GetSymmetryEigenstatesMap()->begin()->second == NULL)
    {
        *errstream << "SymmetryEigenstatesMap was previously generated but Eigenstates have been deleted, regenerate them first." << std::endl;
        exit(1);
    }

    SymmetryEigenstatesMap::const_iterator sem_it;
    for(sem_it = GetSymmetryEigenstatesMap()->begin(); sem_it != GetSymmetryEigenstatesMap()->end(); sem_it++)
    {
        for(int i = 0; i < sem_it->second->GetNumEigenvalues(); i++)
        {
            ConfigGenerator* confgen = GenerateConfigurations(sem_it->first);
            const RelativisticConfigList* configs = confgen->GetRelConfigs();
            RelativisticConfigList::const_iterator list_it = configs->begin();
            std::map<Configuration, double> percentages;
            unsigned int N = 0;
            while(list_it != configs->end())
            {   N += list_it->NumJStates();
                list_it++;
            }

            int j = 0;
            list_it = configs->begin();
            while(list_it != configs->end())
            {
                Configuration nrconfig(list_it->GetNonRelConfiguration());
                if(percentages.find(nrconfig) == percentages.end())
                    percentages[nrconfig] = 0.;
        
                for(unsigned int Jstate = 0; Jstate < list_it->NumJStates(); Jstate++)
                {
                    double coeff = sem_it->second->GetEigenvectors()[i*N + j];
                    coeff = coeff * coeff * 100;
        
                    percentages[nrconfig] += coeff;
                    j++;
                }
        
                list_it++;
            }

            if(sem_it->second->GetgFactors() != NULL)
            {
                GetSolutionMap()->insert(std::pair<SolutionID, Solution>(SolutionID(sem_it->first, i),Solution(sem_it->second->GetEigenvalues()[i], percentages, sem_it->second->GetgFactors()[i])));
            }
            else
            {
                GetSolutionMap()->insert(std::pair<SolutionID, Solution>(SolutionID(sem_it->first, i),Solution(sem_it->second->GetEigenvalues()[i], percentages)));
            }
            delete confgen;
        }
    }
}
