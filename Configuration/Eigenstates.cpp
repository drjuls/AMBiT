#include "Eigenstates.h"
#include "Include.h"

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
