#include "Include.h"
#include "Atom.h"

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
