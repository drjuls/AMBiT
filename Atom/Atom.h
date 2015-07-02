#ifndef ATOM_H
#define ATOM_H

#include "Atom/GetPot"
#include "Configuration/HamiltonianMatrix.h"
#include "HartreeFock/Core.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"
#include "Configuration/Level.h"

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
    /** For all levels requested, create continuum wave and calculate autoionization rate.
        PRE: ChooseHamiltoniansAndRead() must have been run.
     */
    void Autoionization(pLevelConst target);

    /** Use continuum energy grid for autoionisation.
        PRE: ChooseHamiltoniansAndRead() must have been run.
     */
    void AutoionizationEnergyGrid(pLevelConst target);

    /** Configuration averaged autoionisation rate.
        PRE: Have to have nrconfigs.
     */
    void AutoionizationConfigurationAveraged(pLevelConst target);
    void AutoionizationConfigurationAveraged(const OccupationMap& target);

protected:
    /** Initialise AngularDataLibrary if it hasn't been already.
        If trial is included, check that the number of valence electrons is correct before using.
     */
    void InitialiseAngularDataLibrary(pAngularDataLibrary trial = nullptr);

    /** Instantiate LevelMap and read from file.
        Get set of requested symmetries from user input.
     */
    pLevelStore ChooseHamiltonians(pConfigList nrlist);

    /** Get configurations based on single electron configurations with no CI.
        PRE: MakeBasis() must have been run.
     */
    LevelVector SingleElectronConfigurations(pHamiltonianID sym);

    /** Attempt to read basis from file and generate HF operator.
        Return true if successful, false if file "identifier.basis" not found.
     */
    bool ReadBasis();

    void PrintBasis();

public:
    /** Generate integrals with MBPT. CIIntegralsMBPT will automatically store them,
        each processor with a separate file.
        Optionally include delta (cm-1). (eg: delta = E_CI - E_HF for the ground state)
     */
    void InitialiseIntegralsMBPT(bool CoreMBPT = true, bool ValenceMBPT = false);

    /** Collate integrals generated from CIIntegralsMBPT.
        num_processors is the number of stored files to collate.
        If num_processors = 0, it will assume "NumProcessors" (i.e. all processors).
     */
    void CollateIntegralsMBPT(unsigned int num_processors = 0);

public:
    void GenerateCowanInputFile();
    void PrintWavefunctionCowan(FILE* fp, pOrbitalConst ds);
    bool ReadGraspMCDF(const std::string& filename);
    void WriteGraspMCDF() const;
    void WriteGraspMcdfOrbital(FILE* fp, pOrbitalConst ds, unsigned int lattice_size) const;

protected:
    /** Parse basis string definition (e.g. 8spd5f) and convert to vector. */
    bool ParseBasisSize(const char* basis_def, std::vector<unsigned int>& num_states);

    /** Assumes excited object has been created (and excited_mbpt, if UseMBPT is true). */
    void CreateBasis(bool UseMBPT);

protected:
    MultirunOptions user_input;

    std::string identifier;
    double Z;         // Nuclear charge

    unsigned int num_valence_electrons;

    pLattice lattice;
    pCore open_core;                //!< Open-shell HF core
    pHFOperator hf;                 //!< Closed-shell HF operator
    pHFOperator hf_open;            //!< Open-shell HF operator
    pOrbitalManagerConst orbitals;
    pHartreeY hartreeY;

    pOneElectronIntegrals hf_electron;              //!< One-body Hamiltonian
    pTwoElectronCoulombOperator twobody_electron;   //!< Two-body Hamiltonian
    pConfigList nrconfigs;
    pAngularDataLibrary angular_library;
    pLevelStore levels;

//    CIIntegralsMBPT* integralsMBPT;
//    CoreMBPTCalculator* mbpt;
//    ValenceCalculator* valence_mbpt;
//    Sigma3Calculator* sigma3;

//    SymmetryEigenstatesMap symEigenstates;
//    SolutionMap* mSolutionMap;

    // CI + MBPT parameters
    unsigned int NumSolutions;
    double MaxEnergy;
    std::string mbptBasisString;    // User input string
    std::string ciBasisString;      // User input string
    bool check_size_only;           // Check integral and CI sizes
    bool save_eigenstates;          // Save final energy eigenstates
    bool generate_mbpt_integrals;   // Generate MBPT integrals on this run
    // Which MBPT bits to include
    bool includeSigma1;
    bool includeSigma2;
    bool includeSigma3;
};

#endif
