#include "Include.h"
#include "RMatrixPrimer.h"
#include "OutStreams.h"
#include "Debug.h"
#include "Universal/ExpLattice.h"

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
        MPI::Intracomm& comm_world = MPI::COMM_WORLD; // Alias
        NumProcessors = comm_world.Get_size();
        ProcessorRank = comm_world.Get_rank();
    #else
        NumProcessors = 1;
        ProcessorRank = 0;
    #endif

    OutStreams::InitialiseStreams();

    try
    {   RMatrixPrimer A(12, 2, "Mg");
        A.GenerateCore();
        A.CreateBasis();
        A.GetCISpectrum();
        //A.GenerateRMatrixInputFile();
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

void RMatrixPrimer::GenerateCore()
{
    /** Build lattice */

    //lattice = new Lattice(1000, 1.e-6, 50.);
    //lattice = new ExpLattice(257, 1.e-5, 0.05);
    lattice = new ExpLattice(320, 1.e-5, 0.05);

    /** Do relativistic Hartree-Fock. */

    DebugOptions.LogFirstBuild(true);       // Set to false if logging not required
    DebugOptions.LogHFIterations(true);
    //DebugOptions.HartreeEnergyUnits(true);
    core = new Core(lattice, Z, Charge);
    core->Initialise();

    DebugOptions.LogFirstBuild(false);
    DebugOptions.LogHFIterations(false);
}

void RMatrixPrimer::CreateBasis()
{
    /** Create basis for CI */

    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);

    // Choose type of basis. Options include BSplineBasis, RStates, RSinStates, and CustomBasis.
    excited = new BSplineBasis(lattice, core);
    excited->SetIdentifier(&identifier);

    // BSplineBasis needs some extra parameters.
    if(dynamic_cast<BSplineBasis*>(excited))
        dynamic_cast<BSplineBasis*>(excited)->SetParameters(40, 7, 45.);

    // Number of excited states in each partial wave
    std::vector<unsigned int> num_states;
    num_states.push_back(2);
    num_states.push_back(2);
    num_states.push_back(2);
    num_states.push_back(1);

    excited->CreateExcitedStates(num_states);

    DebugOptions.OutputHFExcited(false);
}

void RMatrixPrimer::GetCISpectrum()
{
    // Generate Slater (Coulomb) Integrals
    integrals = new CIIntegralsMBPT(*excited);
    integrals->IncludeValenceSMS(false);
    integrals->SetIdentifier(identifier);
    integrals->Update();

    generator = new ConfigGenerator(excited);

    unsigned int two_j;
    NumSolutions = 6;       // Number of eigenstates to solve for

    *outstream << "\nGS Parity:\n" << std::endl;
    for(two_j = 0; two_j <= 0; two_j += 2)
    {   CalculateEnergy(two_j, even);
    }

    *outstream << "\nOpposite Parity:\n" << std::endl;
    for(two_j = 0; two_j <= 6; two_j+=2)
    {   CalculateEnergy(two_j, odd);
    }
}

void RMatrixPrimer::CalculateEnergy(int twoJ, Parity P)
{
    // Generate set of configurations.
    // Make a set of leading configurations and take promotions from it.
    // Can alternatively make a list of configurations; see ConfigFileGenerator class.
    unsigned int num_promotions = 2;
    std::set<Configuration> leading_configs;

    Configuration config1;
    config1.SetOccupancy(NonRelInfo(3, 0), 2);
    //config1.SetOccupancy(NonRelInfo(4, 0), 1);
    leading_configs.insert(config1);

    //Configuration config2;
    //config2.SetOccupancy(NonRelInfo(3, 0), 1);
    //config2.SetOccupancy(NonRelInfo(3, 1), 1);
    //leading_configs.insert(config2);
    
    generator->ClearConfigLists();
    generator->AddLeadingConfigurations(leading_configs);
    generator->GenerateMultipleExcitationsFromLeadingConfigs(num_promotions, P);
    generator->GenerateRelativisticConfigs();
    generator->GenerateProjections(twoJ);

    HamiltonianMatrix* H;

    #ifdef _MPI
        H = new MPIHamiltonianMatrix(*integrals, generator);
    #else
        H = new HamiltonianMatrix(*integrals, generator);
    #endif

    H->GenerateMatrix();
    H->SolveMatrix(NumSolutions, twoJ, true);
    
    delete H;
}

RMatrixPrimer::RMatrixPrimer(unsigned int atomic_number, int charge, const std::string& atom_identifier):
    Z(atomic_number), Charge(charge), identifier(atom_identifier),
    SD_CI(false), MBPT_CI(false), GenerateFromFile(false), NumSolutions(6),
    excited(NULL), excited_mbpt(NULL), generator(NULL),
    integrals(NULL), integralsMBPT(NULL), mbpt(NULL), sigma3(NULL)
{}

RMatrixPrimer::~RMatrixPrimer(void)
{
    if(generator)
        delete generator;
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

void RMatrixPrimer::GenerateRMatrixInputFile()
{
    // Generate configurations
    Configuration config1;
    config1.SetOccupancy(NonRelInfo(3, 0), 2);
    //leading_configs.insert(config1);
    
    // Get total statistical weight of configurations
    FILE* fp = fopen("Rmatrix.input", "wt");

    // Number of orbitals; points in radial mesh (+1 for R = 0)
    fprintf(fp, " %7d%7d\n ", excited->NumStates(), lattice->Size()+1);
    
    // Print basis states
    int count = 0;
    unsigned int i;
    const DiscreteState* ds;

    // Statistical weights of states
    ConstStateIterator it = excited->GetConstStateIterator();
    it.First();
    while(!it.AtEnd())
    {   if(count == 4)
        {   fprintf(fp, "\n ");
            count = 1;
        }
        else
            count++;
        ds = it.GetState();
        *outstream << ds->Name() << " " << NonRelInfo(ds).Name() << std::endl;
        fprintf(fp, "%16.8E", double(config1.GetOccupancy(NonRelInfo(ds))));
        it.Next();
    }

    // Radial mesh
    fprintf(fp, "\n %16.8E", 0.0);  // Print zero point on lattice
    const double* R = lattice->R();
    count = 1;                      // Already printed one number
    for(i = 0; i < lattice->Size(); i++)
    {   if(count == 4)
        {   fprintf(fp, "\n ");
            count = 1;
        }
        else
            count++;
        fprintf(fp, "%16.8E", R[i]);
    }
    
    it.First();
    while(!it.AtEnd())
    {
        // Upper component
        ds = it.GetState();
        fprintf(fp, "\n %4d  %4d", ds->RequiredPQN(), ds->Kappa());

        fprintf(fp, "\n %16.8E", 0.0);  // Print zero point
        count = 1;
        for(i = 0; i < ds->Size(); i++)
        {   if(count == 4)
            {   fprintf(fp, "\n ");
                count = 1;
            }
            else
                count++;
            fprintf(fp, "%16.8E", ds->f[i]);
        }
        while(i < lattice->Size())
        {   if(count == 4)
            {   fprintf(fp, "\n ");
                count = 1;
            }
            else
                count++;
            fprintf(fp, "%16.8E", 0.0);
            i++;
        }
        
        // Lower component
        ds = it.GetState();
        fprintf(fp, "\n %4d  %4d", ds->RequiredPQN(), ds->Kappa());

        fprintf(fp, "\n %16.8E", 0.0);  // Print zero point
        count = 1;
        for(i = 0; i < ds->Size(); i++)
        {   if(count == 4)
            {   fprintf(fp, "\n ");
                count = 1;
            }
            else
                count++;
            fprintf(fp, "%16.8E", ds->g[i]*Constant::Alpha);
        }
        while(i < lattice->Size())
        {   if(count == 4)
            {   fprintf(fp, "\n ");
                count = 1;
            }
            else
                count++;
            fprintf(fp, "%16.8E", 0.0);
            i++;
        }

        it.Next();
    }
}

void RMatrixPrimer::Write() const
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

void RMatrixPrimer::Read()
{
    // Read core electron states
    std::string filename = identifier + ".core.atom";
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
