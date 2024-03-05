#include "Include.h"
#include "ConfigFileGenerator.h"
#include "Universal/MathConstant.h"

namespace Ambit
{
void ConfigFileGenerator::Clear()
{
    ConfigGenerator::Clear();

    InputContributions.clear();
    Contributions.clear();
}

void ConfigFileGenerator::ReadConfigs(double cutoff)
{
    // Check to see if we haven't restored nrlist since a read.
    if(nrlist.empty() && !rlist.empty())
        RestoreNonRelConfigs();

    FILE* fp = file_err_handler->fopen(input_filename.c_str(), "r");
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
            L = MathConstant::Instance()->GetL(buffer[i]);
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

        if(percentage >= cutoff)
        {   InputContributions[config] = percentage;
            if(config.GetParity() == symmetry.GetParity())
                nrlist.push_back(config);
        }
    }
    
    nrlist.sort();
    nrlist.unique();

    file_err_handler->fclose(fp);
}

void ConfigFileGenerator::WriteConfigs(double cutoff)
{
    std::map<Configuration, double>::const_iterator it = Contributions.begin();
    
    FILE* fp = file_err_handler->fopen(output_filename.c_str(), "w");
    
    while(it != Contributions.end())
    {
        if(it->second >= cutoff)
            fprintf(fp, "%s:\t%5.2lf\n", it->first.Name().c_str(), it->second);
        it++;
    }
    
    file_err_handler->fclose(fp);
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

void ConfigFileGenerator::GenerateMultipleExcitationsFromImportantConfigs(unsigned int num_excitations, double cutoff)
{
    ConfigList important;

    std::set<Configuration>::const_iterator lc_it = leading_configs.begin();    
    while(lc_it != leading_configs.end())
    {
        important.push_back(*lc_it);
        lc_it++;
    }
    
    std::map<Configuration, double>::const_iterator input_it = InputContributions.begin();
    while(input_it != InputContributions.end())
    {
        if(input_it->second >= cutoff)
            important.push_back(input_it->first);
        input_it++;
    }
    
    important.sort();
    important.unique();
    
    GenerateMultipleExcitations(important, num_excitations);
}
}
