#ifndef RMATRIX_PRIMER_H
#define RMATRIX_PRIMER_H

#include "Universal/Constant.h"
#include "HartreeFock/Core.h"
#include "Basis/ExcitedStates.h"
#include "Configuration/Configuration.h"
#include "Configuration/CIIntegrals.h"
#include "Configuration/CIIntegralsMBPT.h"
#include "Configuration/HamiltonianMatrix.h"

class RMatrixPrimer
{
public:
    // Main program
    void CalculateStructure();

public:
    RMatrixPrimer(unsigned int atomic_number, int charge, const std::string& atom_identifier);
    ~RMatrixPrimer();

    /** Identifier is only really used to identify this atom in filenames for storage.
        Should be altered when one of the physical parameters (eg nuclear) changes.
      */
    const std::string& GetIdentifier() { return identifier; }
    void SetIdentifier(const std::string& atom_identifier) { identifier = atom_identifier; }

public:
    /** Generate Hartree-Fock core states. */
    void GenerateCore();

    /** Create a basis of excited states orthogonal to the core states and to each other. */
    void CreateBasis();

    /** Do Configuration Interaction to get energy spectrum. */
    void GetCISpectrum();

    /** Generate a 'target.inp' file for use in R-Matrix code. */
    void GenerateRMatrixInputFile();

public:
    /** Write the atomic state to file, including all known electron wavefunctions. */
    void Write() const;

    /** Read should be called immediately after initialisation of the atom, in order that
        initial Hartree-Fock values are not calculated (just to save computational time).
      */
    void Read();

protected:
    /** Do CI for a particular total ang. momentum and parity. */
    void CalculateEnergy(int twoJ, Parity P); 

private:
    std::string identifier;
    const double Z;         // Nuclear charge
    const double Charge;    // Degree of ionisation

    Lattice* lattice;
    Core* core;
    ExcitedStates* excited;
    ExcitedStates* excited_mbpt;
    ConfigGenerator* generator;
    CIIntegrals* integrals;
    CIIntegralsMBPT* integralsMBPT;
    MBPTCalculator* mbpt;
    Sigma3Calculator* sigma3;

    // Configuration Interaction parameters
    bool SD_CI;     // Only use single and double excitations in CI.
    bool MBPT_CI;   // Include MBPT effects
    bool GenerateFromFile;  // Generate non-relativistic configs from file.
    unsigned int NumSolutions;
};

#endif
