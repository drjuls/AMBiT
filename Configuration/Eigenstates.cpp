#include "Eigenstates.h"
#include "Include.h"
#include "Universal/Constant.h"

Eigenstates::Eigenstates(const std::string& atom_identifier, ConfigGenerator* configlist, bool configlist_owner):
    identifier(atom_identifier), configs(configlist), config_owner(configlist_owner),
    num_eigenvalues(0), num_eigenvectors(0), eigenvalues(NULL), eigenvectors(NULL)
{
    N = configs->GetNumJStates();
}

Eigenstates::Eigenstates(const Eigenstates& other):
    configs(other.configs), N(other.N), config_owner(other.config_owner),
    num_eigenvalues(other.num_eigenvalues), num_eigenvectors(other.num_eigenvectors),
    eigenvalues(other.eigenvalues), eigenvectors(other.eigenvectors)
{}

Eigenstates::~Eigenstates()
{
    if(config_owner)
        delete configs;
    if(num_eigenvalues && eigenvalues)
        delete[] eigenvalues;
    if(num_eigenvectors && eigenvectors)
        delete[] eigenvectors;
}

unsigned int Eigenstates::GetTwoJ() const
{   return configs->GetTwoJ();
}

Parity Eigenstates::GetParity() const
{   return configs->GetParity();
}

ConfigGenerator* Eigenstates::GetConfigGenerator() const
{   return configs;
}

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

void Eigenstates::SetIdentifier(const std::string& atom_identifier)
{   identifier = atom_identifier;
}

void Eigenstates::Write() const
{
    if(ProcessorRank == 0)
    {
        std::string filename = identifier + "." + configs->GetSymmetry().GetString() + ".eigenstates";

        FILE* fp = fopen(filename.c_str(), "wb");

        // Write N as sanity check
        fwrite(&N, sizeof(unsigned int), 1, fp);

        // Write eigenvalues
        fwrite(&num_eigenvalues, sizeof(unsigned int), 1, fp);
        fwrite(eigenvalues, sizeof(double), num_eigenvalues, fp);

        // Write eigenvectors
        fwrite(&num_eigenvectors, sizeof(unsigned int), 1, fp);
        fwrite(eigenvectors, sizeof(double), N * num_eigenvectors, fp);

        fclose(fp);
    }
}

bool Eigenstates::Read()
{
    std::string filename = identifier + "." + configs->GetSymmetry().GetString() + ".eigenstates";

    FILE* fp = fopen(filename.c_str(), "rb");
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

    // Write eigenvectors
    if(eigenvectors)
        delete[] eigenvectors;

    fread(&num_eigenvectors, sizeof(unsigned int), 1, fp);
    eigenvectors = new double[N * num_eigenvectors];
    fread(eigenvectors, sizeof(double), N * num_eigenvectors, fp);

    fclose(fp);
    
    return true;
}

bool Eigenstates::Restore()
{
    bool ret = configs->Read();
    if(ret)
        ret = Read();

    return ret;
}

void Eigenstates::Clear()
{
    if(config_owner)
        configs->Clear();

    if(num_eigenvalues && eigenvalues)
    {   num_eigenvalues = 0;
        delete[] eigenvalues;
    }
    if(num_eigenvectors && eigenvectors)
    {   num_eigenvectors = 0;
        delete[] eigenvectors;
    }
}

void Eigenstates::Print() const
{
    unsigned int i, j;

    *outstream << "Solutions for J = " << double(GetTwoJ())/2. << ", P = ";
    if(GetParity() == even)
        *outstream << "even:" << std::endl;
    else
        *outstream << "odd:" << std::endl;

    const RelativisticConfigList* configlist = configs->GetRelConfigs();

    for(i=0; i<num_eigenvectors; i++)
    {
        unsigned int solution = i;

        *outstream << i << ": " << std::setprecision(8) << eigenvalues[solution] << "    "
            << std::setprecision(12) << eigenvalues[solution]*Constant::HartreeEnergy_cm << " /cm" << std::endl;

        // Get non-rel configuration percentages
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
        
        // Find most important configuration, and print all leading configurations.
        std::map<Configuration, double>::const_iterator it_largest_percentage = percentages.begin();
        double largest_percentage = 0.0;

        std::map<Configuration, double>::const_iterator it = percentages.begin();
        while(it != percentages.end())
        {
            if(it->second > largest_percentage)
            {   it_largest_percentage = it;
                largest_percentage = it->second;
            }

            if(it->second > 1.)
                *outstream << std::setw(20) << it->first.Name() << "  "<< std::setprecision(2)
                    << it->second << "%" << std::endl;
            it++;
        }

        *outstream << std::endl;
    }
}

void Eigenstates::Print(double max_energy) const
{
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
}
