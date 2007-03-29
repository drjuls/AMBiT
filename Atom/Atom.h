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
    void RunOpen();

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
    void CreateRBasis(const StateInfo* ionised = NULL);
    void CreateBSplineBasis(const StateInfo* ionised = NULL);
    void CreateCustomBasis(const StateInfo* ionised = NULL);

    void DoClosedShellSMS(bool include_mbpt = true);
    void DoClosedShellFSModifyR(bool include_mbpt = true);
    void DoClosedShellVolumeShift(bool include_mbpt = true);
    void DoClosedShellAlphaVar(bool include_mbpt = true);

    /** Calculate CI using the existing basis. */
    void OpenShellEnergy(int twoJ, HamiltonianMatrix* H);

    /** Fill out the rlist based on config and twoJ.
        Creates a new Hamiltonian object, which the user must later delete.
     */
    HamiltonianMatrix* CreateHamiltonian(unsigned int twoJ, unsigned int num_particles, Parity P);

    /** Check sizes of matrices before doing full scale calculation. */
    void CheckMatrixSizes();

    /** Calculate specific mass shift. */
    void DoOpenShellSMS(int twoJ, HamiltonianMatrix* H);

    /** Core relaxation term. */
    void SMS_V0(int twoJ, HamiltonianMatrix* H);

    /** Core-valence term. */
    void SMS_V1(int twoJ, HamiltonianMatrix* H);

    /** Valence-valence term (calculates energy and matrix element). */
    void SMS_V2(int twoJ, HamiltonianMatrix* H);

    /** Calculate volume shift. */
    void DoOpenShellFSModifyR(int twoJ, HamiltonianMatrix* H);
    void DoOpenShellVolumeShift(int twoJ, HamiltonianMatrix* H);

    /** Calculate relativistic shift (q-values). */
    void DoOpenShellAlphaVar(int twoJ, HamiltonianMatrix* H);

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
    Eigenstates* GenerateConfigurations(const Symmetry& sym, bool try_read = true);

    /** Calculate energies for all chosen symmetries.
        PRE: Need to have generated all integrals.
     */
    void CalculateEnergies();

public:
    /** Generate multiple integrals with MBPT. */
    void GenerateMulitipleIntegralsMBPT(bool CoreMBPT = true, bool ValenceMBPT = false, double* delta = NULL);

    /** Collate multiple integrals generated from CIIntegralsMBPT. */
    void CollateMultipleIntegralsMBPT(unsigned int num_processors = 0);

    /** Read stored integrals and calculate any remaining without MBPT.
        The MBPT_CI flag tells it whether to use CIIntegralsMBPT or not (CIIntegralsMBPT stores around
        twice as many integrals because of the loss of symmetry due to box diagrams).
     */
    void GenerateMultipleIntegrals(bool MBPT_CI);

    /** Calculate energies for all chosen symmetries. Integrals are generated as needed. */
    void CalculateMultipleEnergies();

public:
    Eigenstates* GetEigenstates(const Symmetry& sym);
    double GetE1MatrixElement(const ElectronInfo& e1, const ElectronInfo& e2) const;
    double GetE1ReducedMatrixElementSquared(const StateInfo& e1, const StateInfo& e2) const;
    ExcitedStates* GetBasis() { return excited; }

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
    MBPTCalculator* mbpt;
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
};

#endif
