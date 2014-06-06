#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Atom.h"
#include "Basis/BasisGenerator.h"

//#include "MBPT/CoreMBPTCalculator.h"
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
    orbitals = nullptr;

    // CI + MBPT parameters
    NumSolutions = user_input("CI/NumSolutions", 6);
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
{}

pCore Atom::MakeBasis(pCoreConst hf_open_core_start)
{
    bool use_read = true;
    if(user_input.search(2, "--clean", "-c"))
        use_read = false;

    bool use_write = true;
    if(user_input.search(2, "--dont-save", "-d"))
        use_write = false;

    if(!use_read || !ReadBasis())
    {
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

        if(user_input.search("HF/--read-grasp0"))
        {   // Read lattice and core and basis orbitals
    //        ReadGraspMCDF("MCDF.DAT");
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
        hf_core = basis_generator.GenerateHFCore(hf_open_core_start);

        // Create Basis
        orbitals = basis_generator.GenerateBasis();

        // Make closed HF operator
        hf = basis_generator.GetHFOperator();
        pCore closed_core(new Core(*orbitals->core));
        hf->SetCore(closed_core);

        // HartreeY operator
        hartreeY = basis_generator.GetHartreeY();

        if(use_write)
        {
            std::string filename = identifier + ".basis";
            orbitals->Write(filename);
        }
    }

    return hf_core;
}

bool Atom::ReadBasis()
{
    // Import lattice and all orbitals
    std::string filename = identifier + ".basis";
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        return false;
    else
        fclose(fp);

    pOrbitalManager modifiable_orbitals(new OrbitalManager(filename));
    lattice = modifiable_orbitals->GetLattice();

    // Generate HF operator
    BasisGenerator basis_generator(lattice, user_input);
    hf = basis_generator.RecreateBasis(modifiable_orbitals);

    orbitals = modifiable_orbitals;
    pCore closed_core(new Core(*orbitals->core));
    hf->SetCore(closed_core);

    hf_core = basis_generator.GetHFCore();

    // HartreeY operator
    hartreeY = basis_generator.GetHartreeY();

    return true;
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
