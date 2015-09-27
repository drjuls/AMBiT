#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Atom.h"
#include "Basis/BasisGenerator.h"
#include "Universal/ExpLattice.h"
#include "MBPT/BruecknerDecorator.h"
#include "HartreeFock/HartreeFocker.h"
#include "HartreeFock/ConfigurationParser.h"

Atom::Atom(const MultirunOptions userInput, unsigned int atomic_number, const std::string& atom_identifier):
    user_input(userInput), Z(atomic_number), identifier(atom_identifier)
{}

Atom::~Atom(void)
{}

pCore Atom::MakeBasis(pCoreConst hf_open_core_start)
{
    bool use_read = true;
    if(user_input.search(2, "--clean", "-c"))
        use_read = false;

    if(!use_read || !ReadBasis())
    {
        int N = user_input("HF/N", -1);  // Number of electrons
        if(N == -1 || N > Z)
        {   *errstream << "Must specify N (number of electrons) with Z >= N" << std::endl;
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
        open_core = basis_generator.GenerateHFCore(hf_open_core_start);
        hf_open = basis_generator.GetOpenHFOperator();

        // Create Basis
        orbitals = basis_generator.GenerateBasis();

        // Make closed HF operator
        hf = basis_generator.GetClosedHFOperator();

        // HartreeY operator
        hartreeY = basis_generator.GetHartreeY();

        std::string filename = identifier + ".basis";
        orbitals->Write(filename);
    }

    return open_core;
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
    hf_open = basis_generator.RecreateBasis(modifiable_orbitals);

    orbitals = modifiable_orbitals;
    hf = basis_generator.GetClosedHFOperator();

    open_core = basis_generator.GetHFCore();

    // HartreeY operator
    hartreeY = basis_generator.GetHartreeY();

    return true;
}

void Atom::GenerateBruecknerOrbitals(bool generate_sigmas)
{
    pBruecknerDecorator brueckner(new BruecknerDecorator(hf));
    bool use_fg = user_input.search("MBPT/Brueckner/--use-lower");
    bool use_gg = user_input.search("MBPT/Brueckner/--use-lower-lower");
    brueckner->IncludeLower(use_fg, use_gg);

    double sigma_start_r = user_input("MBPT/Brueckner/StartPoint", 4.35e-5);
    double sigma_end_r   = user_input("MBPT/Brueckner/EndPoint", 8.0);
    int stride = user_input("MBPT/Brueckner/Stride", 4);
    brueckner->SetMatrixParameters(stride, sigma_start_r, sigma_end_r);

    // Attempt to read all requested kappas
    std::set<int> valence_kappas;
    for(auto& pair: *orbitals->valence)
        valence_kappas.insert(pair.first.Kappa());

    for(int kappa: valence_kappas)
        brueckner->Read(identifier, kappa);

    // Replace hf operator for rest of calculation
    hf = brueckner;

    // Make new sigma potentials if they haven't been read (slowly)
    if(generate_sigmas)
    {
        for(int kappa: valence_kappas)
        {   brueckner->CalculateSigma(kappa, orbitals, hartreeY);
            brueckner->Write(identifier, kappa);
        }
    }

    // And finally change all our valence orbitals to Brueckner orbitals
    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    HartreeFocker hartree_focker(ode_solver);

    // Get scalings or energies
    unsigned int scaling_length = user_input.vector_variable_size("MBPT/Brueckner/Scaling");
    unsigned int energy_scaling_length = user_input.vector_variable_size("MBPT/Brueckner/EnergyScaling");

    // Run this one unless MBPT/Brueckner/EnergyScaling option is used
    if(scaling_length)
    {
        for(int i = 0; i < scaling_length-1; i+=2)
        {
            int kappa = user_input("MBPT/Brueckner/Scaling", 0, i);
            double scale = user_input("MBPT/Brueckner/Scaling", 0.0, i+1);

            brueckner->SetSigmaScaling(kappa, scale);
        }
    }
    else if(energy_scaling_length)
    {   pOrbital brueckner_orbital;

        *outstream << "Brueckner scaling:" << std::endl;

        // Get orbitals to scale and scale them
        for(int i = 0; i < energy_scaling_length-1; i+=2)
        {
            std::string orbital_info = user_input("MBPT/Brueckner/EnergyScaling", "", i);
            OrbitalInfo info(ConfigurationParser::ParseOrbital<OrbitalInfo>(orbital_info));
            double requested_energy = user_input("MBPT/Brueckner/EnergyScaling", 0.0, i+1);

            pOrbital orbital = orbitals->valence->GetState(info);
            if(orbital)
            {
                double initial_energy = orbital->Energy();
                double requested_energy_change = requested_energy - initial_energy;

                // Make a copy of the old orbital
                brueckner_orbital.reset(new Orbital(*orbital));

                double deltaE;
                double scale = 1.;
                do
                {   brueckner->SetSigmaScaling(info.Kappa(), scale);

                    // Iterate
                    hartree_focker.ConvergeOrbitalAndExchange(brueckner_orbital, brueckner, &HartreeFocker::IterateOrbital, 1.e-7 * fabs(requested_energy));

                    // Rescale
                    deltaE = brueckner_orbital->Energy() - initial_energy;
                    scale = requested_energy_change/deltaE * scale;
                } while(!std::isnan(scale) && fabs(scale) < 2. && fabs((deltaE - requested_energy_change)/initial_energy) > 1.e-6);

                *outstream << "  kappa = " << info.Kappa() << ": scaling = " << std::setprecision(6) << scale << std::endl;
            }
            else
                *errstream << "Brueckner: Orbital " << info.Name() << " not found in valence set.";
        }

        *outstream << std::endl;
    }

    pOrbital brueckner_orbital;
    for(auto& pair: *orbitals->valence)
    {
        // Make a copy of the old orbital
        brueckner_orbital.reset(new Orbital(*pair.second));

        // Iterate
        hartree_focker.ConvergeOrbitalAndExchange(brueckner_orbital, brueckner, &HartreeFocker::IterateOrbital, hartree_focker.EnergyTolerance);

        // Copy back to orbital manager
        *pair.second = *brueckner_orbital;
    }

    *outstream << "Brueckner orbitals:\n";
    orbitals->valence->Print();
    *outstream << std::endl;
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
    #ifdef AMBIT_USE_MPI
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
        #ifdef AMBIT_USE_MPI
            MPI::COMM_WORLD.Barrier();
        #endif
        Write();
    }
}
*/
