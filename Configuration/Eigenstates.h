#ifndef EIGENSTATES_H
#define EIGENSTATES_H

#include "ConfigGenerator.h"
#include "RelativisticConfiguration.h"
#include "Symmetry.h"

class Eigenstates
{
    /** Storage class for eigenstates of the Hamiltonian with a particular (J, P) symmetry.
        Eigenstates also takes posession of configlist (the ConfigGenerator) if initialised with
        configlist_owner = true.
     */
public:
    Eigenstates(const std::string& atom_identifier, ConfigGenerator* configlist, bool configlist_owner = true);
    Eigenstates(const Eigenstates& other);
    virtual ~Eigenstates();

    unsigned int GetTwoJ() const;
    Parity GetParity() const;

    ConfigGenerator* GetConfigGenerator() const;

    RelativisticConfigList* GetRelConfigs() const
    {   return configs->GetRelConfigs();
    }

    unsigned int GetNumEigenvalues() const
    {   return num_eigenvalues;
    }

    unsigned int GetNumEigenvectors() const
    {   return num_eigenvectors;
    }

    unsigned int GetEigenvectorLength() const
    {   return N;
    }

    void SetEigenvalues(double* values, unsigned int num_values);
    const double* GetEigenvalues() const;

    void SetEigenvectors(double* vectors, unsigned int num_vectors);
    const double* GetEigenvectors() const;

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
    
protected:
    std::string identifier;

    unsigned int N;             // size of Hamiltonian matrix, length of eigenvectors
    bool config_owner;
    ConfigGenerator* configs;

    unsigned int num_eigenvalues;   // number of eigenvalues stored
    double* eigenvalues;
    unsigned int num_eigenvectors;  // number of eigenvectors stored
    double* eigenvectors;
};

typedef std::map<Symmetry, Eigenstates*> SymmetryEigenstatesMap;

#endif
