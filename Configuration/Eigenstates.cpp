#include "Eigenstates.h"
#include "Include.h"
#include "Universal/MathConstant.h"

#define CONFIG_PRINT_LIMIT 5.

namespace Ambit
{
Eigenstates::Eigenstates(const std::string& atom_identifier, pRelativisticConfigListConst configlist):
    identifier(atom_identifier), configs(configlist),
    num_eigenvalues(0), num_eigenvectors(0), num_gFactors(0), eigenvalues(nullptr), eigenvectors(nullptr), gFactors(nullptr)
{
    N = configs->NumCSFs();
}

Eigenstates::Eigenstates(const Eigenstates& other):
    configs(other.configs), N(other.N), config_owner(other.config_owner),
    num_eigenvalues(other.num_eigenvalues), num_eigenvectors(other.num_eigenvectors), num_gFactors(other.num_gFactors),
    eigenvalues(other.eigenvalues), eigenvectors(other.eigenvectors), gFactors(other.gFactors)
{}

Eigenstates::~Eigenstates()
{
    if(num_eigenvalues && eigenvalues)
        delete[] eigenvalues;
    if(num_eigenvectors && eigenvectors)
        delete[] eigenvectors;
    if(num_gFactors && gFactors)
        delete[] gFactors;
}

//unsigned int Eigenstates::GetTwoJ() const
//{   return configs->GetTwoJ();
//}
//
//Parity Eigenstates::GetParity() const
//{   return configs->GetParity();
//}

void Eigenstates::SetEigenvalues(double* values, unsigned int num_values)
{
    if(eigenvalues)
        delete[] eigenvalues;

    num_eigenvalues = num_values;
    eigenvalues = values;
}

const double* Eigenstates::GetEigenvalues() const
{   return eigenvalues;
}

void Eigenstates::SetEigenvectors(double* vectors, unsigned int num_vectors)
{
    if(eigenvectors)
        delete[] eigenvectors;

    num_eigenvectors = num_vectors;
    eigenvectors = vectors;
}

const double* Eigenstates::GetEigenvectors() const
{   return eigenvectors;
}

void Eigenstates::SetgFactors(double* g_factors, unsigned int num_gfactors)
{
    if(gFactors)
        delete[] gFactors;
    num_gFactors = num_gfactors;
    gFactors = g_factors;
}

const double* Eigenstates::GetgFactors() const
{   return gFactors;
}

void Eigenstates::SetIdentifier(const std::string& atom_identifier)
{   identifier = atom_identifier;
}

void Eigenstates::Write() const
{
    if(ProcessorRank == 0)
    {
//        std::string filename = identifier + "." + configs->GetSymmetry().GetString() + ".eigenstates";
        std::string filename = identifier + ".eigenstates";

        FILE* fp = file_err_handler->fopen(filename.c_str(), "wb");

        // Write N as sanity check
        file_err_handler->fwrite(&N, sizeof(unsigned int), 1, fp);

        // Write eigenvalues
        file_err_handler->fwrite(&num_eigenvalues, sizeof(unsigned int), 1, fp);
        file_err_handler->fwrite(eigenvalues, sizeof(double), num_eigenvalues, fp);

        // Write eigenvectors
        file_err_handler->fwrite(&num_eigenvectors, sizeof(unsigned int), 1, fp);
        file_err_handler->fwrite(eigenvectors, sizeof(double), N * num_eigenvectors, fp);

        // Write g-factors
        if(num_gFactors && gFactors)
        {   file_err_handler->fwrite(&num_gFactors, sizeof(unsigned int), 1, fp);
            file_err_handler->fwrite(gFactors, sizeof(double), num_gFactors, fp);
        }
        else
        {   unsigned int zero = 0;
            file_err_handler->fwrite(&zero, sizeof(unsigned int), 1, fp);
        }

        file_err_handler->fclose(fp);
    }
}

bool Eigenstates::Read()
{
//    std::string filename = identifier + "." + configs->GetSymmetry().GetString() + ".eigenstates";
    std::string filename = identifier + ".eigenstates";

    FILE* fp = file_err_handler->fopen(filename.c_str(), "rb");
    if(!fp)
        return false;

    // Read N and check for consistency
    unsigned int size;
    fread(&size, sizeof(unsigned int), 1, fp);
    if(N != size)
    {   *errstream << "Eigenstates::Read: Stored number of JStates doesn't match:\n"
                   << "   N = " << N << ", stored size = " << size << std::endl;
        return false;
    }

    // Read eigenvalues
    if(eigenvalues)
        delete[] eigenvalues;

    fread(&num_eigenvalues, sizeof(unsigned int), 1, fp);
    eigenvalues = new double[num_eigenvalues];
    fread(eigenvalues, sizeof(double), num_eigenvalues, fp);

    // Read eigenvectors
    if(eigenvectors)
        delete[] eigenvectors;

    fread(&num_eigenvectors, sizeof(unsigned int), 1, fp);
    eigenvectors = new double[N * num_eigenvectors];
    fread(eigenvectors, sizeof(double), N * num_eigenvectors, fp);

    // Read g-factors
    if(gFactors)
        delete[] gFactors;

    fread(&num_gFactors, sizeof(unsigned int), 1, fp);
    if(num_gFactors)
    {   gFactors = new double[num_gFactors];
        fread(gFactors, sizeof(double), num_gFactors, fp);
    }
    else
        gFactors = nullptr;

    file_err_handler->fclose(fp);
    
    return true;
}

bool Eigenstates::Restore()
{
    // Does not need to be restored
    if(num_eigenvalues && num_eigenvectors)
        return true;

    bool ret = configs->Read();
    if(ret)
        ret = Read();

    return ret;
}

void Eigenstates::Clear()
{
    if(num_eigenvalues && eigenvalues)
    {   num_eigenvalues = 0;
        delete[] eigenvalues;
    }
    if(num_eigenvectors && eigenvectors)
    {   num_eigenvectors = 0;
        delete[] eigenvectors;
    }
    if(gFactors)
        delete[] gFactors;
}

void Eigenstates::Print() const
{
    *outstream << "Solutions for J = " << double(GetTwoJ())/2. << ", P = ";
    if(GetParity() == even)
        *outstream << "even:" << std::endl;
    else
        *outstream << "odd:" << std::endl;

    unsigned int i = 0;
    while(i < num_eigenvectors)
    {
        Print(i, CONFIG_PRINT_LIMIT);
        *outstream << std::endl;
        i++;
    }
}

void Eigenstates::Print(double max_energy) const
{
    *outstream << "Solutions for J = " << double(GetTwoJ())/2. << ", P = ";
    if(GetParity() == even)
        *outstream << "even:" << std::endl;
    else
        *outstream << "odd:" << std::endl;

    unsigned int i = 0;
    while((i < num_eigenvectors) && (eigenvalues[i] <= max_energy))
    {
        Print(i, CONFIG_PRINT_LIMIT);
        *outstream << std::endl;
        i++;
    }
}

void Eigenstates::Print(unsigned int solution, double config_fraction_limit) const
{
    if(solution >= num_eigenvalues)
        return;

    unsigned int j;
    const RelativisticConfigList* configlist = configs->GetRelConfigs();

    *outstream << solution << ": " << std::setprecision(8) << eigenvalues[solution] << "    "
        << std::setprecision(12) << eigenvalues[solution]*MathConstant::Instance()->HartreeEnergyInInvCm() << " /cm" << std::endl;

    // Get non-rel configuration percentages
    if(solution < num_eigenvectors)
    {
        RelativisticConfigList::const_iterator list_it = configlist->begin();
        std::map<Configuration, double> percentages;  // Map non-rel configurations to percentages

        j = 0;
        while(list_it != configlist->end())
        {
            Configuration nrconfig(list_it->GetNonRelConfiguration());
            if(percentages.find(nrconfig) == percentages.end())
                percentages[nrconfig] = 0.;

            for(unsigned int Jstate = 0; Jstate < list_it->NumJStates(); Jstate++)
            {
                double coeff = eigenvectors[solution*N + j];
                coeff = coeff * coeff * 100;

                percentages[nrconfig] += coeff;
                j++;
            }

            list_it++;
        }
        
        // Print leading configurations.
        if(config_fraction_limit)
        {
            std::map<Configuration, double>::const_iterator it = percentages.begin();
            while(it != percentages.end())
            {
                if(it->second > config_fraction_limit)
                    *outstream << std::setw(20) << it->first.Name() << "  "<< std::setprecision(2)
                        << it->second << "%" << std::endl;
                it++;
            }
        }
        // Print only major contribution.
        else
        {
            std::map<Configuration, double>::const_iterator it_largest_percentage = percentages.begin();
            std::map<Configuration, double>::const_iterator it_second_largest_percentage = percentages.begin();
            double largest_percentage = 0.0;
            double second_largest_percentage = 0.0;

            std::map<Configuration, double>::const_iterator it = percentages.begin();
            while(it != percentages.end())
            {
                if(it->second > largest_percentage)
                {   // largest -> second_largest
                    it_second_largest_percentage = it_largest_percentage;
                    second_largest_percentage = largest_percentage;
                    // it -> largest
                    it_largest_percentage = it;
                    largest_percentage = it->second;
                }
                else if(it->second > second_largest_percentage)
                {   // it -> second_largest
                    it_second_largest_percentage = it;
                    second_largest_percentage = it->second;
                }

                it++;
            }

            if(!second_largest_percentage || (largest_percentage/second_largest_percentage > 2.))
            {   // Leading config only
                *outstream << std::setw(20) << it_largest_percentage->first.Name() << "  " << std::setprecision(2)
                    << largest_percentage << "%" << std::endl;
            }
            else if(largest_percentage + second_largest_percentage > 80.)
            {   // Leading two configs only
                *outstream << std::setw(20) << it_largest_percentage->first.Name() << "  "
                    << std::setprecision(2) << largest_percentage << "%" << std::endl;
                *outstream << std::setw(20) << it_second_largest_percentage->first.Name() << "  "
                    << std::setprecision(2) << second_largest_percentage << "%" << std::endl;
            }
            else
            {   // Print all configurations with more than 5% contribution
                std::map<Configuration, double>::const_iterator it = percentages.begin();
                while(it != percentages.end())
                {
                    if(it->second > CONFIG_PRINT_LIMIT)
                        *outstream << std::setw(20) << it->first.Name() << "  "<< std::setprecision(2)
                            << it->second << "%" << std::endl;
                    it++;
                }
            }
        }

        if(gFactors)
            *outstream << "    g-factor = " << std::setprecision(5) << gFactors[solution] << std::endl;
    }
}

void Eigenstates::PrintCowan(FILE* fp, double energy_shift) const
{
    unsigned int solution;
    double E_average = 0.0;
    for(solution=0; solution<num_eigenvectors; solution++)
    {
        E_average += eigenvalues[solution];
    }
    E_average = E_average/num_eigenvectors;

    fprintf(fp, "%5d%14.7E\n", num_eigenvectors, E_average + energy_shift);

    unsigned int i;
    for(solution=0; solution<num_eigenvectors; solution++)
    {
        int P;
        if(configs->GetParity() == even)
            P = 1;
        else
            P = -1;

        fprintf(fp, "%14.7E %4d %4d\n", eigenvalues[solution] + energy_shift, configs->GetTwoJ()+1, P);

        unsigned int count = 0;

        for(i = 0; i < N; i++)
        {   if(count == 5)
            {   fprintf(fp, "\n");
                count = 1;
            }
            else
                count++;
            fprintf(fp, "%14.7E", eigenvectors[solution*N + i]);
        }
        fprintf(fp, "\n");
    }

    const RelativisticConfigList* configlist = configs->GetRelConfigs();
    RelativisticConfigList::const_iterator it = configlist->begin();

    fprintf(fp, "\nConfigurations:\n");
    i = 1;
    while(it != configlist->end())
    {   fprintf(fp, "%4d %s\n", i, it->Name().c_str());
        it++; i++;
    }
    fprintf(fp, "\n\n");
}

double SymmetryEigenstatesMap::GetLargestTwoJ(Parity aParity)
{
    SymmetryEigenstatesMap::iterator it = begin();
    double largest = 0;
    while(it != end())
    {
        if(it->first.GetParity() == aParity && it->first.GetTwoJ() > largest)
        {
            largest = it->first.GetTwoJ();
        }
        it++;
    }
    return largest;
}

bool SymmetryEigenstatesMap::RestoreAll()
{
    bool result = true;
    SymmetryEigenstatesMap::iterator it = begin();
    while(it != end())
    {
        if(it->second != NULL)
        {
            result = result && it->second->Restore();
        }
        it++;
    }
    
    return result;
}
}
