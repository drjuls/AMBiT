#include "Include.h"
#include "ConfigFileGenerator.h"
#include "Universal/Constant.h"
    
void ConfigFileGenerator::ReadConfigs(Parity parity, double cutoff)
{
    FILE* fp = fopen(input_filename.c_str(), "r");
    char buffer[200];

    int L, pqn, occupancy;
    char tempbuf[20];
    unsigned int i, j;

    while(fgets(buffer, 200, fp))
    {
        i = 0;
        Configuration config;
        double percentage;

        // Move to first configuration
        while(!isdigit(buffer[i]) && buffer[i] != '\0')
        {   i++;
        }

        while(buffer[i] != '\0')
        {
            j=0;

            // Get config
            while(isdigit(buffer[i]))
            {   tempbuf[j] = buffer[i];
                i++; j++;
            }
            tempbuf[j] = 0;
            pqn = atoi(tempbuf);

            // Get L
            for(L = 0; L < 10; L++)
            {   if(Constant::SpectroscopicNotation[L] == buffer[i])
                    break;
            }
            i++;

            // Get occupancy
            j = 0;
            while(isdigit(buffer[i]))
            {   tempbuf[j] = buffer[i];
                i++; j++;
            }
            tempbuf[j] = 0;
            occupancy = atoi(tempbuf);

            config.SetOccupancy(NonRelInfo(pqn, L), occupancy);
            
            // End of core?
            if(buffer[i] == ':')
            {   i++;
                break;
            }

            // Next configuration term
            while(!isdigit(buffer[i]) && buffer[i] != '\0')
                i++;
        }
        
        percentage = atof(buffer+i);

        if((percentage >= cutoff) && (config.GetParity() == parity))
        {   nrlist.push_back(config);
        }
    }
    
    nrlist.sort();
    nrlist.unique();

    fclose(fp);
}

void ConfigFileGenerator::WriteConfigs(double cutoff)
{
    std::map<Configuration, double>::const_iterator it = Contributions.begin();
    
    FILE* fp = fopen(output_filename.c_str(), "w");
    
    while(it != Contributions.end())
    {
        if(it->second >= cutoff)
            fprintf(fp, "%s:\t%5.2lf\n", it->first.Name().c_str(), it->second);
        it++;
    }
    
    fclose(fp);
}

void ConfigFileGenerator::AddPercentages(const std::map<Configuration, double> percentages)
{
    std::map<Configuration, double>::const_iterator it = percentages.begin();
    
    while(it != percentages.end())
    {
        if(Contributions.find(it->first) == Contributions.end())
            Contributions[it->first] = it->second;
        else
            Contributions[it->first] = mmax(Contributions[it->first], it->second);

        it++;
    }
}


