#include "CustomBasis.h"
#include "Include.h"

void CustomBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    FILE* fp = fopen("CustomBasis.txt", "r");
    if(fp == NULL)
    {   printf("Unable to open file CustomBasis.txt");
        PAUSE;
        exit(1);
    }

    char buffer[100];

    char tempbuf[20];
    unsigned int i, j;

    while(fgets(buffer, 100, fp))
    {
        // Get state details
        i = 0;
        NonRelInfo final = ReadNonRelInfo(buffer, i);

        // Get method of state creation
        j = 0;
        while(!isspace(buffer[i]))
        {   tempbuf[j] = buffer[i];
            i++; j++;
        }
        tempbuf[j] = 0;

        if(!strcmp(tempbuf, "HF"))
        {
            DiscreteState* ds = new DiscreteState(lattice, final.PQN(), final.GetFirstRelativisticInfo().Kappa());
            unsigned int it = core->CalculateExcitedState(ds);
            if(it)
                Orthogonalise(ds);
            AddState(ds);

            unsigned int count, r;
            double fmax = 0.;
            for(count = 0; count < ds->Size(); count++)
                if(fabs(ds->f[count]) > fmax)
                {   fmax = ds->f[count];
                    r = count;
                }
            if(core->GetDebugOptions().DebugHFExcited())
                printf("  %s  en: %f  size: %d  Rmax: %d  NZ: %d\n", ds->Name().c_str(), ds->Energy(), ds->Size(),
                            r, int(ds->RequiredPQN()) - int(ds->NumZeroes()) - int(ds->L()) - 1);
            //fmax = 0.;
            //CoupledFunction exch;
            //core->CalculateExchange(*ds, exch);
            //for(count = 0; count < ds->Size(); count++)
            //    fmax += fabs(ds->df[count]/lattice->dR(count) + double(ds->Kappa())/lattice->R(count)*ds->f[count]
            //    - (2. + Constant::AlphaSquared * core->GetHFPotential()[count])*ds->g[count]
            //    - Constant::AlphaSquared*exch.g[count]);
            //printf("    Deviation: %f\n", fmax);

            if(final.L() != 0)
            {   ds = new DiscreteState(lattice, final.PQN(), final.GetSecondRelativisticInfo().Kappa());
                unsigned int it = core->CalculateExcitedState(ds);
                if(it)
                    Orthogonalise(ds);
                AddState(ds);
                fmax = 0.;
                for(count = 0; count < ds->Size(); count++)
                    if(fabs(ds->f[count]) > fmax)
                    {   fmax = ds->f[count];
                        r = count;
                    }
                if(core->GetDebugOptions().DebugHFExcited())
                    printf("  %s  en: %f  size: %d  Rmax: %d  NZ: %d\n", ds->Name().c_str(), ds->Energy(), ds->Size(),
                            r, int(ds->RequiredPQN()) - int(ds->NumZeroes()) - int(ds->L()) - 1);
            }
        }
        else
        {   NonRelInfo prev = ReadNonRelInfo(buffer, i);
            if((!strcmp(tempbuf, "RS") && (prev.L() + 1 != final.L())) ||
               (strcmp(tempbuf, "RS") && (prev.L() != final.L())))
            {   fprintf(stderr, "Error in \"CustomBasis.txt\", %s -> %s.\n", prev.Name().c_str(), final.Name().c_str());
                PAUSE
                exit(1);
            }
            DiscreteState* ds = new DiscreteState(lattice, final.PQN(), final.GetFirstRelativisticInfo().Kappa());
            DiscreteState* previous = GetState(prev.GetFirstRelativisticInfo().GetStateInfo());
            if(previous == NULL)
                previous = core->GetState(prev.GetFirstRelativisticInfo().GetStateInfo());
            if(previous == NULL)
            {   printf("Error in \"CustomBasis.txt\": %s undefined.\n", prev.Name().c_str());
                PAUSE
                exit(1);
            }
            if(!strcmp(tempbuf, "R"))
                MultiplyByR(previous, ds);
            else if(!strcmp(tempbuf, "S"))
                MultiplyBySinR(previous, ds);
            else if(!strcmp(tempbuf, "RS"))
                MultiplyByRSinR(previous, ds);
            else
            {   printf("Error in \"CustomBasis.txt\": %s unknown.\n", tempbuf);
                PAUSE
                exit(1);
            }
            AddState(ds);
            unsigned int count, r;
            double fmax = 0.;
            for(count = 0; count < ds->Size(); count++)
                if(fabs(ds->f[count]) > fmax)
                {   fmax = ds->f[count];
                    r = count;
                }
            if(core->GetDebugOptions().DebugHFExcited())
                printf("  %s  en: %f  size: %d  Rmax: %d  NZ: %d\n", ds->Name().c_str(), ds->Energy(), ds->Size(),
                            r, int(ds->RequiredPQN()) - int(ds->NumZeroes()) - int(ds->L()) - 1);
            //fmax = 0.;
            //for(count = 0; count < ds->Size(); count++)
            //    fmax += fabs(ds->df[count]/lattice->dR(count) + ds->Kappa()/lattice->R(count)*ds->f[count]
            //    -(2. + Constant::AlphaSquared * (core->GetHFPotential()[count] + core->GetLocalExchangeApproximation()[count]))*ds->g[count]);
            //printf("    Deviation: %f\n", fmax);
            
            if(final.L() != 0)
            {   ds = new DiscreteState(lattice, final.PQN(), final.GetSecondRelativisticInfo().Kappa());
                previous = GetState(prev.GetSecondRelativisticInfo().GetStateInfo());
                if(previous == NULL)
                    previous = core->GetState(prev.GetSecondRelativisticInfo().GetStateInfo());
                if(previous == NULL)
                {   printf("Error in \"CustomBasis.txt\": %s undefined.\n", prev.Name().c_str());
                    PAUSE
                    exit(1);
                }
                if(!strcmp(tempbuf, "R"))
                    MultiplyByR(previous, ds);
                else if(!strcmp(tempbuf, "S"))
                    MultiplyBySinR(previous, ds);
                else if(!strcmp(tempbuf, "RS"))
                    MultiplyByRSinR(previous, ds);
                AddState(ds);
                fmax = 0.;
                for(count = 0; count < ds->Size(); count++)
                    if(fabs(ds->f[count]) > fmax)
                    {   fmax = ds->f[count];
                        r = count;
                    }
                if(core->GetDebugOptions().DebugHFExcited())
                    printf("  %s  en: %f  size: %d  Rmax: %d  NZ: %d\n", ds->Name().c_str(), ds->Energy(), ds->Size(),
                            r, int(ds->RequiredPQN()) - int(ds->NumZeroes()) - int(ds->L()) - 1);
                //if(core->GetDebugOptions().DebugHFExcited())
                //    std::cout << "  " << ds->Name() << " en:   " << ds->Energy() << "  size:  " << ds->Size() << std::endl;
            }
        }
    }

    fclose(fp);
}

/** Update all of the excited states because the core has changed. */
void CustomBasis::Update()
{
    Clear();
    CreateExcitedStates(std::vector<unsigned int>());
}

NonRelInfo CustomBasis::ReadNonRelInfo(char* buffer, unsigned int& num_read)
{
    while(isspace(buffer[num_read]))
        num_read++;

    char tempbuf[10];
    unsigned int j = 0;
    while(isdigit(buffer[num_read]))
    {   tempbuf[j] = buffer[num_read];
        num_read++; j++;
    }
    tempbuf[j] = 0;
    int pqn = atoi(tempbuf);

    // Get L
    unsigned int L;
    for(L = 0; L < 10; L++)
    {   if(Constant::SpectroscopicNotation[L] == buffer[num_read])
            break;
    }
    num_read++;

    while(isspace(buffer[num_read]) && (buffer[num_read] != '\n'))
        num_read++;

    return NonRelInfo(pqn, L);
}