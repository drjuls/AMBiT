#ifndef ATOM_H
#define ATOM_H

#include "Universal/Constant.h"
#include "HartreeFock/Core.h"
#include "Basis/ExcitedStates.h"
#include "Configuration/Configuration.h"
#include "Configuration/CIIntegrals.h"
#include "Configuration/HamiltonianMatrix.h"

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

    /** Create a second order sigma, if it doesn't already exist. */
    void GetSigma(const StateInfo& info);

    void DoClosedShellSMS(bool include_mbpt = true);
    void DoClosedShellVolumeShift(bool include_mbpt = true);
    void DoClosedShellAlphaVar(bool include_mbpt = true);

    /** Calculate CI using the existing basis.
        If size_only, then do not calculate CI, but report the size of the matrix.
     */
    void OpenShellEnergy(int twoJ, const Configuration& config, bool size_only = false);
    void OpenShellEnergy(int twoJ, HamiltonianMatrix* H);

    /** Fill out the rlist based on config and twoJ.
        Creates a new Hamiltonian object, which the user must later delete.
     */
    HamiltonianMatrix* CreateHamiltonian(int twoJ, const Configuration& config, RelativisticConfigList& rlist);

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
    void DoOpenShellVolumeShift(int twoJ, HamiltonianMatrix* H);

    /** Calculate relativistic shift (q-values). */
    void DoOpenShellAlphaVar(int twoJ, HamiltonianMatrix* H);

private:
    std::string identifier;
    const double Z;         // Nuclear charge
    const double Charge;    // Degree of ionisation

    Lattice* lattice;
    Core* core;
    ExcitedStates* excited;
    CIIntegrals* integrals;

    // Configuration Interaction parameters
    bool SD_CI;     // Only use single and double excitations in CI.
    unsigned int NumSolutions;
};

#endif
