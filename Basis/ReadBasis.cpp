#include "Include.h"
#include "ReadBasis.h"
#include "Universal/Constant.h"
#include "Universal/Interpolator.h"
#include "HartreeFock/StateIntegrator.h"
#include "Universal/ExpLattice.h"
#include <cctype>

void ReadBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    FILE* fp = fopen(filename.c_str(), "r");
    if(fp == NULL)
    {   *errstream << "Unable to open file " << filename << std::endl;
        PAUSE;
        exit(1);
    }

    // Lattice of stored states
    Lattice* read_lattice = new ExpLattice(300, 1.e-5, 0.05);

    std::vector<unsigned int> num_states_so_far(num_states_per_l.size(), 0);

    unsigned int i;
    unsigned int total_num_states = 0;
    unsigned int num_states_read = 0;

    for(i = 0; i < num_states_per_l.size(); i++)
    {   total_num_states += num_states_per_l[i];
        if(i > 0)
            total_num_states += num_states_per_l[i];
    }

    StateIntegrator I(lattice);

    while(num_states_read < total_num_states)
    {
        // Read upper component
        int pqn, l, kappa;
        char lchar;
        int numpoints;
        double energy;

        fscanf(fp, "%d", &pqn);
        fscanf(fp, "%c", &lchar);
        lchar = tolower(lchar);
        for(l = 0; l < 10; l++)
            if(Constant::SpectroscopicNotation[l] == lchar)
                break;
        if(l >= 10)
        {   *errstream << "unable to find l from notation " << lchar << std::endl;
            PAUSE;
            exit(1);
        }
        
        lchar = getc(fp);
        if(lchar == '-')
            kappa = l;
        else
            kappa = -(l+1);
            
        fscanf(fp, "%d", &numpoints);
        fscanf(fp, "%*12le%*12le");     // Coefficients of r^Kappa
        fscanf(fp, "%12le", &energy);   // Energy relative to threshold (eV)

        DiscreteState* ds;
        DiscreteState* ds_readlattice = new DiscreteState(pqn, kappa);
        ds_readlattice->ReSize(numpoints);

        for(i = 0; i<numpoints; i++)
            fscanf(fp, "%14le", &(ds_readlattice->f[i]));

        *outstream << pqn << " " << l << lchar << "." << numpoints
                   << " " << energy << std::endl;

        *outstream << ds_readlattice->f[0] << " " << ds_readlattice->f[numpoints-1] << std::endl;

        // Read lower component
        fscanf(fp, "%d", &pqn);
        fscanf(fp, "%c", &lchar);
        lchar = tolower(lchar);
        for(l = 0; l < 10; l++)
            if(Constant::SpectroscopicNotation[l] == lchar)
                break;
        lchar = getc(fp);
        if(lchar == '-')
            kappa = l;
        else
            kappa = -(l+1);
        fscanf(fp, "%d", &numpoints);
        fscanf(fp, "%*12le%*12le");

        if(pqn != ds_readlattice->RequiredPQN() || kappa != ds_readlattice->Kappa() ||
           numpoints != ds_readlattice->Size())
        {   *errstream << "Lower component doesn't match upper in " << ds_readlattice->Name() << std::endl;
            PAUSE;
            exit(1);
        }

        double g;
        for(i = 0; i<numpoints; i++)
        {   fscanf(fp, "%14le", &g);
            ds_readlattice->g[i] = g/Constant::Alpha;
        }

        *outstream << ds_readlattice->g[0]*Constant::Alpha << " "
                   << ds_readlattice->g[numpoints-1]*Constant::Alpha << std::endl;

        unsigned int order = 6;
        if(*read_lattice == *lattice)
        {   // Get derivative
            ds = ds_readlattice;
            Interpolator interp(lattice);

            interp.GetDerivative(ds->f, ds->df, order);
            interp.GetDerivative(ds->g, ds->dg, order);
            *outstream << "No interp." << std::endl;
        }
        else
        {   // Interpolate onto regular lattice
            Interpolator interp(read_lattice);            
            ds = new DiscreteState(pqn, kappa);

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
            *outstream << "Interp." << std::endl;
        }

        AddState(ds);
        num_states_read++;

        ds->SetEnergy(I.HamiltonianMatrixElement(*ds, *ds, *core));
    }

    fclose(fp);

    if(DebugOptions.OutputHFExcited())
        *outstream << "Basis Orthogonality test: " << TestOrthogonality() << std::endl;
}
