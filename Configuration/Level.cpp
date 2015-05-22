#include "Level.h"
#include "Include.h"
#include "Universal/MathConstant.h"
#include "NonRelConfiguration.h"

HamiltonianID::HamiltonianID(const std::string& name):
    sym(-1)
{
    int twoJ;
    std::stringstream ss(name);
    ss >> twoJ;

    Parity P;
    char c_parity;
    ss.get(c_parity);
    if(c_parity == 'o')
        P = Parity::odd;
    else
        P = Parity::even;

    sym = Symmetry(twoJ, P);
}

std::string HamiltonianID::Name() const
{
    std::string ret = itoa(sym.GetTwoJ());
    ret.append(LetterName(sym.GetParity()));

    return ret;
}

std::string HamiltonianID::Print() const
{
    std::string ret = "J = " + boost::lexical_cast<std::string>(sym.GetJ());
    ret += ", P = " + LowerName(sym.GetParity());

    return ret;
}

void HamiltonianID::Write(FILE* fp) const
{
    int two_j = sym.GetTwoJ();
    Parity P = sym.GetParity();
    fwrite(&two_j, sizeof(int), 1, fp);
    fwrite(&P, sizeof(Parity), 1, fp);
    
    // Config list
    configs->Write(fp);
}

void HamiltonianID::Read(FILE* fp)
{
    int two_j;
    Parity P;
    fread(&two_j, sizeof(int), 1, fp);
    fread(&P, sizeof(Parity), 1, fp);
    sym = Symmetry(two_j, P);

    // Config list
    pRelativisticConfigList new_configs(std::make_shared<RelativisticConfigList>());
    new_configs->Read(fp);
    configs = new_configs;
}

Level::Level(const double& energy, const double* csf_eigenvector, pHamiltonianID hamiltonian_id, unsigned int numCSFs):
    eigenvalue(energy), eigenvector(numCSFs), hamiltonian(hamiltonian_id)
{
    unsigned int N = numCSFs;
    if(numCSFs == 0)
    {   N = hamiltonian->GetRelativisticConfigList()->NumCSFs();
        eigenvector.resize(N);
    }

    memcpy(eigenvector.data(), csf_eigenvector, N * sizeof(double));
    gFactor = std::nan("");
}

const Level& Level::operator=(const Level& other)
{
    hamiltonian = other.hamiltonian;
    eigenvalue = other.eigenvalue;
    eigenvector = other.eigenvector;
    gFactor = other.gFactor;

    return *this;
}

Level& Level::operator=(Level&& other)
{
    hamiltonian.swap(other.hamiltonian);
    eigenvalue = other.eigenvalue;
    gFactor = other.gFactor;
    eigenvector.swap(other.eigenvector);

    return *this;
}
