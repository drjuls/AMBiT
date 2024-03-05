#ifndef ATOM_H
#define ATOM_H

#include "GetPot"
#include "Configuration/HamiltonianMatrix.h"
#include "HartreeFock/Core.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"
#include "Configuration/LevelMap.h"
#include "HartreeFock/NucleusDecorator.h"
#include "MBPT/Sigma3Calculator.h"

namespace Ambit
{
class Sigma3Calculator;

/** Atom is something like a "LevelGenerator", but with Read/Write operations,
    multiprocessor functions, etc, based on user input.
 */
class Atom
{
public:
    Atom(const MultirunOptions userInput, unsigned int atomic_number, const std::string& atom_identifier);
    ~Atom();

    /** Read existing basis or perform Hartree-Fock calculation
        (optionally using hf_open_core_start as a starting approximation).
        Set up basis orbitals and operators; write orbitals.
        Return open-shell core orbitals.
     */
    pCore MakeBasis(pCoreConst hf_open_core_start = nullptr);

    /** Get the open-shell core orbitals. */
    pCore GetOpenShellCore() { return open_core; }

    /** Get the orbital basis. */
    pOrbitalManagerConst GetBasis() const { return orbitals; }

    /** Get the lattice. */
    pLattice GetLattice() { return lattice; }

    /** Get physical constants.
        PRE: MakeBasis() must have been run.
     */
    pPhysicalConstant GetPhysicalConstants() { return hf->GetPhysicalConstant(); }
    pPhysicalConstantConst GetPhysicalConstants() const { return hf->GetPhysicalConstant(); }

    /** Get closed-shell HF operator. */
    pHFOperatorConst GetHFOperator() const { return hf; }
    pHFOperator GetHFOperator() { return hf; }

    /** Get open-shell HF operator. */
    pHFOperatorConst GetOpenHFOperator() const { return hf_open; }
    pHFOperator GetOpenHFOperator() { return hf_open; }

    /** Get HartreeY operator. */
    pHartreeY GetHartreeY() { return hartreeY; }

    /** Get nuclear if it exists. */
    pNucleusDecorator GetNucleusDecorator() { return nucleus; }

public:
    /** Generate integrals with MBPT, store and collate from all processors.
        PRE: MakeBasis() must have been run.
     */
    void MakeMBPTIntegrals();

    /** Generate integrals for CI using any stored MBPT integrals and calculating the rest with no MBPT.
        PRE: MakeBasis() must have been run.
     */
    void MakeIntegrals();

    /** Clear integrals for CI. */
    void ClearIntegrals();

    /** Check sizes of matrices before doing full scale calculation.
        PRE: MakeBasis() must have been run.
     */
    void CheckMatrixSizes(pAngularDataLibrary angular_lib = nullptr);

    /** Get levels to be calculated, so the client can request one be done at a time.
        PRE: MakeBasis() must have been run.
     */
    pLevelStore ChooseHamiltoniansAndRead(pAngularDataLibrary angular_lib = nullptr);

    /** Calculate levels for all chosen symmetries in user input.
        PRE: ChooseHamiltoniansAndRead() must have been run.
     */
    pLevelStore CalculateEnergies();

    /** Calculate levels for given HamiltonianID, generating CI integrals if required via. MakeIntegrals().
        PRE: ChooseHamiltoniansAndRead() must have been run.
        Levels are added to LevelMap levels.
    */
    LevelVector CalculateEnergies(pHamiltonianID hID);

    pLevelStore GetLevels() { return levels; }
    pAngularDataLibrary GetAngularDataLibrary() { return angular_library; }

public:
    /** For all levels requested, create continuum wave and calculate autoionization width.
        PRE: ChooseHamiltoniansAndRead() must have been run.
     */
    void Autoionization(const LevelVector& target);

    /** Use continuum energy grid for autoionisation.
        PRE: ChooseHamiltoniansAndRead() must have been run.
     */
    void AutoionizationEnergyGrid(const LevelVector& target);

    /** Configuration averaged autoionisation strength S.
        PRE: Have to have nrconfigs.
     */
    void AutoionizationConfigurationAveraged(const LevelVector& target);
    void AutoionizationConfigurationAveraged(const OccupationMap& target);

    void InternalConversion(const LevelVector& target);
    void InternalConversionConfigurationAveraged(const LevelVector& target);

protected:
    /** Initialise AngularDataLibrary if it hasn't been already.
        If trial is included, check that the number of valence electrons is correct before using.
     */
    void InitialiseAngularDataLibrary(pAngularDataLibrary trial = nullptr);

    /** Instantiate LevelMap and read from file.
        Get set of requested symmetries from user input.
     */
    pLevelStore ChooseHamiltonians(pRelativisticConfigList rlist);

    /** Get configurations based on single electron configurations with no CI.
        PRE: MakeBasis() must have been run.
     */
    LevelVector SingleElectronConfigurations(pHamiltonianID sym);

    /** Attempt to read basis from file and generate HF operator.
        Return true if successful, false if file "identifier.basis" not found.
     */
    bool ReadBasis();

    /** Generate or read sigma matrices, create Brueckner orbitals, and use them everywhere.
        If(!generate_sigmas), new sigmas will not be created.
     */
    void GenerateBruecknerOrbitals(bool generate_sigmas);

public:
    void GenerateCowanInputFile();
    void PrintWavefunctionCowan(FILE* fp, pOrbitalConst ds);
    bool ReadGraspMCDF(const std::string& filename);
    void WriteGraspMCDF() const;
    void WriteGraspMcdfOrbital(FILE* fp, pOrbitalConst ds, unsigned int lattice_size) const;

protected:
    MultirunOptions user_input;

    std::string identifier;
    double Z;         // Nuclear charge

    pLattice lattice;
    pCore open_core;                //!< Open-shell HF core
    pHFOperator hf;                 //!< Closed-shell HF operator
    pHFOperator hf_open;            //!< Open-shell HF operator
    pOrbitalManagerConst orbitals;
    pHartreeY hartreeY;
    pNucleusDecorator nucleus;      //!< Pointer to nucleus (if included)

    pHFIntegrals hf_electron;                       //!< One-body Hamiltonian operator
    pTwoElectronCoulombOperator twobody_electron;   //!< Two-body Hamiltonian operator
    pSigma3Calculator threebody_electron;           //!< Three-body Hamiltonian operator

    pConfigList leading_configs;
    pRelativisticConfigList allconfigs;
    pAngularDataLibrary angular_library;
    pLevelStore levels;
};

}
#endif
