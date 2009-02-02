#ifndef ATOM_H
#define ATOM_H

#include "Universal/Constant.h"
#include "HartreeFock/Core.h"
#include "Basis/ExcitedStates.h"
#include "Configuration/Configuration.h"
#include "Configuration/CIIntegrals.h"
#include "Configuration/CIIntegralsMBPT.h"
#include "Configuration/HamiltonianMatrix.h"
#include "Configuration/Eigenstates.h"

class Atom
{
public:
    // Main program
    void Run();
    void RunOpen(double radius = 0.0);

public:
    Atom(unsigned int atomic_number, int charge, const std::string& atom_identifier, bool read = false);
    ~Atom();

    /** Identifier is only really used to identify this atom in filenames for storage.
        Should be altered when one of the physical parameters (eg nuclear) changes.
      */
    const std::string& GetIdentifier() { return identifier; }
    void SetIdentifier(const std::string& atom_identifier) { identifier = atom_identifier; }

    /** Write the atomic state to file, including all known electron wavefunctions. */
    void Write() const;

    /** Read should be called immediately after initialisation of the atom, in order that
        initial Hartree-Fock values are not calculated (just to save computational time).
      */
    void Read();

public:
    /** Get energy of state given kappa and principal quantum number. */
    double GetEnergy(const StateInfo& info);

    /** TODO: Create a "complete" Hartree-Fock basis, including discrete and continuum states. */

    /** Create a virtual basis, which is discrete yet takes into account parts of the continuum. */
    void CreateRBasis(bool UseMBPT = false);
    void CreateBSplineBasis(bool UseMBPT = false);
    void CreateBSplineBasis(double radius);
    void CreateCustomBasis(bool UseMBPT = false);
    void CreateHartreeFockBasis(bool UseMBPT = false);

    void DoClosedShellSMS(bool include_mbpt = true);
    void DoClosedShellFSModifyR(bool include_mbpt = true);
    void DoClosedShellVolumeShift(bool include_mbpt = true);
    void DoClosedShellAlphaVar(bool include_mbpt = true);

public:
    /** Generate integrals with MBPT. CIIntegralsMBPT will automatically store them,
        each processor with a separate file.
        Optionally include delta (cm-1). (eg: delta = E_CI - E_HF for the ground state)
     */
    void GenerateIntegralsMBPT(bool CoreMBPT = true, bool ValenceMBPT = false, double delta = 0.0);

    /** Collate integrals generated from CIIntegralsMBPT.
        num_processors is the number of stored files to collate.
        If num_processors = 0, it will assume "NumProcessors" (i.e. all processors).
     */
    void CollateIntegralsMBPT(unsigned int num_processors = 0);

    /** Read stored integrals and calculate any remaining without MBPT. */
    void GenerateIntegrals(bool MBPT_CI);

    /** Fill symmetries in Symmetry-Eigenstates map.*/
    void ChooseSymmetries();
    
    /** Generate configurations, projections, and JStates for a symmetry, and store on disk.
        if(try_read), then it checks to see if the configurations exist on file already.
     */
    ConfigGenerator* GenerateConfigurations(const Symmetry& sym, bool try_read = true);

    /** Check sizes of matrices before doing full scale calculation. */
    void CheckMatrixSizes();

    /** Calculate energies for all chosen symmetries.
        PRE: Need to have generated all integrals.
     */
    void CalculateEnergies();

public:
    /** Run atomic energy code multiple times, varying some parameter like
        alpha, NuclearInverseMass (for SMS), or volume shift parameter.
     */
    void RunMultiple(bool include_mbpt = true, bool closed_shell = false);

    /** Set up core, excited states, and integrals some value of the varying parameter.
        PRE: Need to have initialised CIIntegrals* integrals, core and excited states.
     */
    void SetMultipleIntegralsAndCore(unsigned int index);

    /** Generate multiple integrals with MBPT. */
    void GenerateMultipleIntegralsMBPT(bool CoreMBPT = true, bool ValenceMBPT = false, double* delta = NULL);

    /** Collate multiple integrals generated from CIIntegralsMBPT. */
    void CollateMultipleIntegralsMBPT(unsigned int num_processors = 0);

    /** Read stored integrals and calculate any remaining without MBPT.
        The MBPT_CI flag tells it whether to use CIIntegralsMBPT or not (CIIntegralsMBPT stores around
        twice as many integrals because of the loss of symmetry due to box diagrams).
     */
    void GenerateMultipleIntegrals(bool MBPT_CI);

    /** Calculate energies for all chosen symmetries. Integrals are generated as needed. */
    void CalculateMultipleEnergies();

    void CalculateMultipleClosedShell(bool include_mbpt);

public:
    /** Get previously calculated Eigenstates of Symmetry sym. */
    Eigenstates* GetEigenstates(const Symmetry& sym);

    /** Get the set of atomic core orbitals. */
    Core* GetCore() { return core; }

    /** Get the set of valence orbitals. */
    ExcitedStates* GetBasis() { return excited; }

public:
    void GenerateCowanInputFile();
    void PrintWavefunctionCowan(FILE* fp, const DiscreteState* ds);

protected:
    /** Assumes excited object has been created (and excited_mbpt, if UseMBPT is true). */
    void CreateBasis(bool UseMBPT);

private:
    std::string identifier;
    const double Z;         // Nuclear charge
    const double Charge;    // Degree of ionisation

    Lattice* lattice;
    Core* core;
    ExcitedStates* excited;
    ExcitedStates* excited_mbpt;
    CIIntegrals* integrals;
    CIIntegralsMBPT* integralsMBPT;
    CoreMBPTCalculator* mbpt;
    ValenceCalculator* valence_mbpt;
    Sigma3Calculator* sigma3;

    SymmetryEigenstatesMap SymEigenstates;

    // Configuration Interaction parameters
    bool SD_CI;     // Only use single and double excitations in CI.
    bool MBPT_CI;   // Include MBPT effects
    unsigned int NumSolutions;

    bool multiple_SMS;
    bool multiple_alpha;
    bool multiple_volume;
    bool multiple_radius;
    std::vector<double> multiple_parameters;
    std::vector<std::string> multiple_ids;
    std::vector<double> mbpt_delta;

    double alpha0;  // Original value of alpha.
};

#endif
