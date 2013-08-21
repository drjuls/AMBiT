#include "Include.h"
#include "ReadBasis.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/Interpolator.h"
#include "HartreeFock/StateIntegrator.h"
#include "Universal/ExpLattice.h"
#include <cctype>

// Control different types of headers, etc.

// Grasp0 beginning
//  - total number of orbitals    lattice size
//  - occupation numbers for all orbitals
//  - lattice points
#define GRASP0_FILE_START  1

// grasp0 type orbital header:
//      pqn     kappa
#define GRASP0_FILE  1

// cmccore type orbital header:
//      nl[-]   numpoints   coeffs. of r^Kappa     E(eV)
//  eg:   2P-       277   0.0000E+00  0.0000E+00  4.5027E+00
#define CMCCORE_FILE  0

void ReadBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    CreateExcitedStates(num_states_per_l, NULL);
}

void ReadBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l, Core* atom_core)
{
    FILE* fp = fopen(filename.c_str(), "r");
    if(fp == NULL)
    {   *errstream << "Unable to open file " << filename << std::endl;
        PAUSE;
        exit(1);
    }

    Clear();

    // Lattice of stored states
    Lattice* read_lattice = new ExpLattice(303, 1.e-5, 0.05);

    // Need to treat all kappas
    std::vector<unsigned int> num_states_so_far(2 * num_states_per_l.size() - 1, 0);

    unsigned int i;
    double buf;
    unsigned int total_num_states = 0;
    unsigned int num_states_read = 0;
    int numpoints;

    if(GRASP0_FILE_START)
    {
        fscanf(fp, "%d", &total_num_states);
        fscanf(fp, "%d", &numpoints);
        
        // Read occupation numbers
        double total_num_electrons = 0;
        for(i = 0; i < total_num_states; i++)
        {   fscanf(fp, "%le", &buf);
            total_num_electrons += buf;
        }
        if(fabs(total_num_electrons - core->GetZ() + core->GetCharge()) > 1.e-4)
        {   *errstream << "ReadBasis: Sum of occupation numbers not equal to number of electrons." << std::endl;
            *errstream << "   Sum = " << total_num_electrons << std::endl;
            exit(1);
        }

        // Read lattice points
        for(i = 0; i < numpoints; i++)
            fscanf(fp, "%le", &buf);

        // Skip initial zero of all wavefunctions
        numpoints--;
    }
    else
    {   for(i = 0; i < num_states_per_l.size(); i++)
        {   total_num_states += num_states_per_l[i];
            if(i > 0)
                total_num_states += num_states_per_l[i];
        }
    }

    StateIntegrator I(lattice);

    while(num_states_read < total_num_states)
    {
        // Read upper component
        int pqn, l, kappa;
        char lchar;
        double energy;

        if(CMCCORE_FILE)
        {
            fscanf(fp, "%d", &pqn);
            fscanf(fp, "%c", &lchar);
            lchar = tolower(lchar);
            l = MathConstant::Instance()->GetL(lchar);
            
            lchar = getc(fp);
            if(lchar == '-')
                kappa = l;
            else
                kappa = -(l+1);
            
            fscanf(fp, "%d", &numpoints);
            fscanf(fp, "%*12le%*12le");     // Coefficients of r^Kappa
            fscanf(fp, "%12le", &energy);   // Energy relative to threshold (eV)
        }
        else
        {   fscanf(fp, "%d", &pqn);
            fscanf(fp, "%d", &kappa);
        }

        Orbital* ds;
        Orbital* ds_readlattice = new Orbital(pqn, kappa);
        ds_readlattice->ReSize(numpoints);

        if(CMCCORE_FILE)
        {   for(i = 0; i<numpoints; i++)
                fscanf(fp, "%14le", &(ds_readlattice->f[i]));
        }
        else
        {   // skip initial 0.0
            fscanf(fp, "%le", &buf);
            for(i = 0; i<numpoints; i++)
                fscanf(fp, "%le", &(ds_readlattice->f[i]));
        }

        *logstream << ds_readlattice->Name() << ": f read, " << std::flush;

        // Read lower component
        if(CMCCORE_FILE)
        {
            fscanf(fp, "%d", &pqn);
            fscanf(fp, "%c", &lchar);
            lchar = tolower(lchar);
            l = MathConstant::Instance()->GetL(lchar);
            lchar = getc(fp);
            if(lchar == '-')
                kappa = l;
            else
                kappa = -(l+1);
            fscanf(fp, "%d", &numpoints);
            fscanf(fp, "%*12le%*12le");
        }
        else
        {   fscanf(fp, "%d", &pqn);
            fscanf(fp, "%d", &kappa);
        }

        if(pqn != ds_readlattice->RequiredPQN() || kappa != ds_readlattice->Kappa() ||
           numpoints != ds_readlattice->Size())
        {   *errstream << "Lower component doesn't match upper in " << ds_readlattice->Name() << std::endl;
            PAUSE;
            exit(1);
        }

        double g;
        if(CMCCORE_FILE)
        {   for(i = 0; i<numpoints; i++)
            {   fscanf(fp, "%14le", &g);
                ds_readlattice->g[i] = g;
            }
        }
        else
        {   // skip initial 0.0
            fscanf(fp, "%le", &buf);
            for(i = 0; i<numpoints; i++)
            {   fscanf(fp, "%le", &g);
                ds_readlattice->g[i] = g;
            }
        }

        *logstream << "g read" << std::endl;

        unsigned int order = 6;
        if(*read_lattice == *lattice)
        {   // Get derivative
            ds = ds_readlattice;
            Interpolator interp(lattice);

            interp.GetDerivative(ds->f, ds->df, order);
            interp.GetDerivative(ds->g, ds->dg, order);
        }
        else
        {   // Interpolate onto regular lattice
            Interpolator interp(read_lattice);            
            ds = new Orbital(pqn, kappa);

            unsigned int size = lattice->real_to_lattice(read_lattice->R(ds_readlattice->Size()-1));
            ds->ReSize(size);

            const double* R = lattice->R();
            const double* dR = lattice->dR();

            double dfdr, dgdr;
            for(unsigned int i = 0; i < ds->Size(); i++)
            {
                interp.Interpolate(ds_readlattice->f, R[i], ds->f[i], dfdr, order);
                interp.Interpolate(ds_readlattice->g, R[i], ds->g[i], dgdr, order);
                ds->df[i] = dfdr * dR[i];
                ds->dg[i] = dgdr * dR[i];
            }
            
            delete ds_readlattice;
        }

        // Check whether this is a core state
        const Orbital* existing = core->GetState(OrbitalInfo(pqn, kappa));

        // If state is not in the open shell part, check whether it is in the core
        if(existing != NULL && !core->IsOpenShellState(OrbitalInfo(pqn, kappa)))
        {
            // Replace core states if atom_core != NULL
            if(atom_core)
            {   Orbital* ds_replace = atom_core->GetState(OrbitalInfo(pqn, kappa));
                *ds_replace = *ds;
            }
        }
        else
        {   unsigned int index;
            if(kappa > 0)
                index = 2*kappa - 1;
            else
                index = 2*abs(kappa) - 2;

            if(num_states_so_far[index] < num_states_per_l[ds->L()])
                AddState(ds);

            num_states_so_far[index]++;
        }

        num_states_read++;

        ds->SetEnergy(I.HamiltonianMatrixElement(*ds, *ds, *core));
    }

    fclose(fp);

    if(DebugOptions.OutputHFExcited())
        *outstream << "Basis Orthogonality test: " << TestOrthogonality() << std::endl;
}
