#include "CustomBasis.h"
#include "Include.h"

void CustomBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    FILE* fp = fopen("CustomBasis.txt", "r");
    if(fp == NULL)
    {   *errstream << "Unable to open file CustomBasis.txt" << std::endl;
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
            if(DebugOptions.OutputHFExcited())
                *outstream << "  " << ds->Name() << "  en: " << std::setprecision(8) << ds->Energy()
                           << "  size: " << ds->Size() << "  Rmax: " << r << "  NZ: " 
                           << int(ds->RequiredPQN()) - int(ds->NumZeroes()) - int(ds->L()) - 1 << std::endl;

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
                if(DebugOptions.OutputHFExcited())
                    *outstream << "  " << ds->Name() << "  en: " << std::setprecision(8) << ds->Energy()
                               << "  size: " << ds->Size() << "  Rmax: " << r << "  NZ: " 
                               << int(ds->RequiredPQN()) - int(ds->NumZeroes()) - int(ds->L()) - 1 << std::endl;
            }
        }
        else
        {   NonRelInfo prev = ReadNonRelInfo(buffer, i);
            if((!strcmp(tempbuf, "RS") && (prev.L() + 1 != final.L())) ||
               (strcmp(tempbuf, "RS") && (prev.L() != final.L())))
            {   *errstream << "Error in \"CustomBasis.txt\", " << prev.Name() << " -> " << final.Name() << std::endl;
                PAUSE
                exit(1);
            }
            DiscreteState* ds = new DiscreteState(lattice, final.PQN(), final.GetFirstRelativisticInfo().Kappa());
            const DiscreteState* previous = GetState(prev.GetFirstRelativisticInfo().GetStateInfo());
            if(previous == NULL)
                previous = core->GetState(prev.GetFirstRelativisticInfo().GetStateInfo());
            if(previous == NULL)
            {   *errstream << "Error in \"CustomBasis.txt\": " << prev.Name() << " undefined." << std::endl;
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
            {   *errstream << "Error in \"CustomBasis.txt\": " << tempbuf << " unknown." << std::endl;
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
            if(DebugOptions.OutputHFExcited())
                *outstream << "  " << ds->Name() << "  en: " << std::setprecision(8) << ds->Energy()
                           << "  size: " << ds->Size() << "  Rmax: " << r << "  NZ: " 
                           << int(ds->RequiredPQN()) - int(ds->NumZeroes()) - int(ds->L()) - 1 << std::endl;
            
            if(final.L() != 0)
            {   ds = new DiscreteState(lattice, final.PQN(), final.GetSecondRelativisticInfo().Kappa());
                previous = GetState(prev.GetSecondRelativisticInfo().GetStateInfo());
                if(previous == NULL)
                    previous = core->GetState(prev.GetSecondRelativisticInfo().GetStateInfo());
                if(previous == NULL)
                {   *outstream << "Error in \"CustomBasis.txt\": " << prev.Name() << " undefined." << std::endl;
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
                if(DebugOptions.OutputHFExcited())
                    *outstream << "  " << ds->Name() << "  en: " << std::setprecision(8) << ds->Energy()
                               << "  size: " << ds->Size() << "  Rmax: " << r << "  NZ: " 
                               << int(ds->RequiredPQN()) - int(ds->NumZeroes()) - int(ds->L()) - 1 << std::endl;
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
