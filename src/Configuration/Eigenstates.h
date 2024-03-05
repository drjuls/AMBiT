#ifndef EIGENSTATES_H
#define EIGENSTATES_H

#include <map>
#include "ConfigGenerator.h"
#include "RelativisticConfiguration.h"
#include "Symmetry.h"

namespace Ambit
{
/** Storage class for eigenstates of the Hamiltonian with a particular (J, P) symmetry. */
class Eigenstates
{
public:
    Eigenstates(const std::string& atom_identifier, pRelativisticConfigListConst configlist);
    Eigenstates(const Eigenstates& other);
    virtual ~Eigenstates();

//    unsigned int GetTwoJ() const;
//    Parity GetParity() const;

    pRelativisticConfigListConst GetRelConfigs() const
    {   return configs;
    }

    unsigned int GetNumEigenvalues() const
    {   return num_eigenvalues;
    }

    unsigned int GetNumEigenvectors() const
    {   return num_eigenvectors;
    }

    unsigned int GetNumGFactors() const
    {   return num_gFactors;
    }

    unsigned int GetEigenvectorLength() const
    {   return N;
    }

    void SetEigenvalues(double* values, unsigned int num_values);
    const double* GetEigenvalues() const;

    void SetEigenvectors(double* vectors, unsigned int num_vectors);
    const double* GetEigenvectors() const;

    void SetgFactors(double* g_factors, unsigned int num_gfactors);
    const double* GetgFactors() const;

    void SetIdentifier(const std::string& atom_identifier);

    /** Store the eigenvalues and eigenvectors.
        Filename is "identifier.twoJ.P.eigenstates".
     */
    void Write() const;

    /** Read eigenvalues and eigenvectors. Return success. */
    bool Read();

    /** Read configs, eigenvalues, and eigenvectors. Return success. */
    bool Restore();

    /** Clear eigenvalues, and eigenvectors, freeing the memory.
        If config_owner, clear configs too.
     */
    void Clear();

    void Print() const;
    void Print(double max_energy) const;
    void Print(unsigned int solution, double config_fraction_limit = 0.0) const;
    
    void PrintCowan(FILE* fp, double energy_shift = 0.0) const;
protected:
    std::string identifier;

    unsigned int N;             // size of Hamiltonian matrix, length of eigenvectors
    bool config_owner;
    pRelativisticConfigListConst configs;

    unsigned int num_eigenvalues;   // number of eigenvalues stored
    double* eigenvalues;
    unsigned int num_eigenvectors;  // number of eigenvectors stored
    double* eigenvectors;
    unsigned int num_gFactors;
    double* gFactors;
};

class SymmetryEigenstatesMap : public std::map<Symmetry, Eigenstates*>
{
public:
    double GetLargestTwoJ(Parity aParity);
    bool RestoreAll();
};

}
#endif
