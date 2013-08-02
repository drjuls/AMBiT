#include "Include.h"
#include "CustomBasis.h"
#include "Universal/MathConstant.h"

void CustomBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    FILE* fp = fopen(filename.c_str(), "r");
    if(fp == NULL)
    {   *errstream << "Unable to open file " << filename << std::endl;
        PAUSE;
        exit(1);
    }

    NumStatesPerL = num_states_per_l;
    std::vector<unsigned int> num_states_so_far(num_states_per_l.size(), 0);

    char buffer[100];

    char tempbuf[20];
    unsigned int i, j;

    while(fgets(buffer, 100, fp))
    {
        // Get state details
        i = 0;
        NonRelInfo final = ReadNonRelInfo(buffer, i);
        
        // Check whether this state is needed.
        // A "null" num_states_per_l means just do all of them.
        if(!NumStatesPerL.size() ||
            ( (final.L() < NumStatesPerL.size()) &&
              (num_states_so_far[final.L()] < NumStatesPerL[final.L()])))
        {
            if(NumStatesPerL.size())
                num_states_so_far[final.L()]++;

            // Get method of state creation
            j = 0;
            while(!isspace(buffer[i]))
            {   tempbuf[j] = buffer[i];
                i++; j++;
            }
            tempbuf[j] = 0;

            if(!strcmp(tempbuf, "HF"))
            {
                if(!core->IsOpenShellState(final.GetFirstRelativisticInfo()) && (core->GetState(final.GetFirstRelativisticInfo()) != NULL))
                {   *errstream << "Error in " << filename << ": " << final.Name() << " is in HF core." << std::endl;
                    exit(1);
                }

                Orbital* ds = GetState(final.GetFirstRelativisticInfo());
                if(ds == NULL)
                {   ds = new Orbital(final.PQN(), final.GetFirstRelativisticInfo().Kappa());
                    unsigned int loops = core->CalculateExcitedState(ds);
                    if(loops)
                        Orthogonalise(ds);
                    AddState(ds);
                }
                else
                {   unsigned int loops = core->UpdateExcitedState(ds);
                    if(loops)
                        Orthogonalise(ds);
                }

                if(DebugOptions.OutputHFExcited())
                {   unsigned int count, r;
                    double fmax = 0.;
                    for(count = 0; count < ds->Size(); count++)
                        if(fabs(ds->f[count]) > fmax)
                        {   fmax = fabs(ds->f[count]);
                            r = count;
                        }

                    *outstream << "  " << ds->Name() << "  en: " << std::setprecision(8) << ds->Energy()
                               << "  size: " << ds->Size() << "  Rmax: " << r << "  NZ: " 
                               << int(ds->RequiredPQN()) - int(ds->NumNodes()) - int(ds->L()) - 1 << std::endl;
                }

                if(final.L() != 0)
                {
                    ds = GetState(final.GetSecondRelativisticInfo());
                    if(ds == NULL)
                    {   ds = new Orbital(final.PQN(), final.GetSecondRelativisticInfo().Kappa());
                        unsigned int loops = core->CalculateExcitedState(ds);
                        if(loops)
                            Orthogonalise(ds);
                        AddState(ds);
                    }
                    else
                    {   unsigned int loops = core->UpdateExcitedState(ds);
                        if(loops)
                            Orthogonalise(ds);
                    }

                    if(DebugOptions.OutputHFExcited())
                    {   unsigned int count, r;
                        double fmax = 0.;
                        for(count = 0; count < ds->Size(); count++)
                            if(fabs(ds->f[count]) > fmax)
                            {   fmax = fabs(ds->f[count]);
                                r = count;
                            }

                        *outstream << "  " << ds->Name() << "  en: " << std::setprecision(8) << ds->Energy()
                                   << "  size: " << ds->Size() << "  Rmax: " << r << "  NZ: " 
                                   << int(ds->RequiredPQN()) - int(ds->NumNodes()) - int(ds->L()) - 1 << std::endl;
                    }
                }
            }
            else    // Not HF
            {
                NonRelInfo prev = ReadNonRelInfo(buffer, i);
//                if((!strcmp(tempbuf, "RS") && (prev.L() + 1 != final.L())) ||
//                   (strcmp(tempbuf, "RS") && (prev.L() != final.L())))
//                {   *errstream << "Error in " << filename << ", " << prev.Name() << " -> " << final.Name() << std::endl;
//                    PAUSE
//                    exit(1);
//                }

                Orbital* ds = GetState(final.GetFirstRelativisticInfo());
                if(ds == NULL)
                {   ds = new Orbital(final.PQN(), final.GetFirstRelativisticInfo().Kappa());
                    AddState(ds);
                }

                const Orbital* previous = GetState(prev.GetFirstRelativisticInfo());
                if(previous == NULL)
                    previous = core->GetState(prev.GetFirstRelativisticInfo());
                if(previous == NULL)
                {   *errstream << "Error in " << filename << ": " << prev.Name() << " undefined." << std::endl;
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
                {   *errstream << "Error in " << filename << ": " << tempbuf << " unknown." << std::endl;
                    PAUSE
                    exit(1);
                }

                if(DebugOptions.OutputHFExcited())
                {   unsigned int count, r;
                    double fmax = 0.;
                    for(count = 0; count < ds->Size(); count++)
                        if(fabs(ds->f[count]) > fmax)
                        {   fmax = fabs(ds->f[count]);
                            r = count;
                        }

                    *outstream << "  " << ds->Name() << "  en: " << std::setprecision(8) << ds->Energy()
                               << "  size: " << ds->Size() << "  Rmax: " << r << "  NZ: " 
                               << int(ds->RequiredPQN()) - int(ds->NumNodes()) - int(ds->L()) - 1 << std::endl;
                }

                if(final.L() != 0)
                {
                    ds = GetState(final.GetSecondRelativisticInfo());
                    if(ds == NULL)
                    {    ds = new Orbital(final.PQN(), final.GetSecondRelativisticInfo().Kappa());
                         AddState(ds);
                    }

                    previous = GetState(prev.GetSecondRelativisticInfo());
                    if(previous == NULL)
                        previous = core->GetState(prev.GetSecondRelativisticInfo());
                    if(previous == NULL)
                    {   *outstream << "Error in " << filename << ": " << prev.Name() << " undefined." << std::endl;
                        PAUSE
                        exit(1);
                    }

                    if(!strcmp(tempbuf, "R"))
                        MultiplyByR(previous, ds);
                    else if(!strcmp(tempbuf, "S"))
                        MultiplyBySinR(previous, ds);
                    else if(!strcmp(tempbuf, "RS"))
                        MultiplyByRSinR(previous, ds);

                    if(DebugOptions.OutputHFExcited())
                    {   unsigned int count, r;
                        double fmax = 0.;
                        for(count = 0; count < ds->Size(); count++)
                            if(fabs(ds->f[count]) > fmax)
                            {   fmax = fabs(ds->f[count]);
                                r = count;
                            }

                        *outstream << "  " << ds->Name() << "  en: " << std::setprecision(8) << ds->Energy()
                                   << "  size: " << ds->Size() << "  Rmax: " << r << "  NZ: " 
                                   << int(ds->RequiredPQN()) - int(ds->NumNodes()) - int(ds->L()) - 1 << std::endl;
                    }
                }
            }
        }
    }

    // Repeat orthogonalisation
    StateIterator it = GetStateIterator();
    it.First();
    while(!it.AtEnd())
    {   Orthogonalise(it.GetState());
        it.Next();
    }

    fclose(fp);
    if(DebugOptions.OutputHFExcited())
        *outstream << "Basis Orthogonality test: " << TestOrthogonality() << std::endl;
}

void CustomBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l, const std::string& file)
{
    SetFile(file);
    CreateExcitedStates(num_states_per_l);
}

/** Update all of the excited states because the core has changed. */
void CustomBasis::Update()
{
    CreateExcitedStates(NumStatesPerL);
}

void CustomBasis::SetFile(const std::string& file)
{
    filename = file;
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
    unsigned int L = MathConstant::Instance()->GetL(buffer[num_read]);
    num_read++;

    while(isspace(buffer[num_read]) && (buffer[num_read] != '\n'))
        num_read++;

    return NonRelInfo(pqn, L);
}
