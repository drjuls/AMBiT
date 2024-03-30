#include "Configuration/HamiltonianMatrix.h"
#include "gtest/gtest.h"
#include "Include.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Basis/BasisGenerator.h"
#include "Configuration/ConfigGenerator.h"
#include "Atom/MultirunOptions.h"
#include "MBPT/OneElectronMBPT.h"
#include "MBPT/CoreValenceIntegrals.h"

using namespace Ambit;

/* Comparison with MBPT corrections in Beloy et al, 2008.
   This tests the core-valence MBPT (no CI) for a relatively large core.*/
TEST(CoreValenceIntegralsTester, CsGroundState)
{
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);

    pLattice lattice(new Lattice(5000, 1.e-6, 50.));

    // CsI - we want to set up the MBPT options, even though we'll only use them for one component
    std::string user_input_string = std::string() + 
        "NuclearRadius=5.675\n" +
        "Z = 55\n" +
        "-s12\n" +
        "[Lattice]\n" +
        "NumPoints=5000\n" +
        "[HF]\n" +
        "N = 54\n" +
        "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6:'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis=6s\n" +
        "[Basis/BSpline]\n" +
        "K=10\n" +
        "[CI]\n" +
        "LeadingConfigurations = '6s1'\n" +
        "ElectronExcitations = 0\n" +
        "HoleExcitations = 0\n" +
        "NumSolutions = 1\n" +
        "EvenParityTwoJ = '1'\n" +
        "[MBPT]\n" +
        "Basis=30spdfghi\n";
    
    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");
    std::string identifier = "Cs_gtest";

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();

    // Generate MBPT integrals
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pCoulombOperator coulomb(new CoulombOperator(lattice));
    pHartreeY hartreeY(new HartreeY(hf->GetIntegrator(), coulomb));
    // Bare integrals for MBPT
    pHFIntegrals bare_one_body_integrals = std::make_shared<HFIntegrals>(orbitals, hf);
    pSlaterIntegrals bare_two_body_integrals = std::make_shared<SlaterIntegralsFlatHash>(orbitals, hartreeY);

    // MBPT calculators
    std::string fermi_orbitals = "";
    pCoreMBPTCalculator core_mbpt = std::make_shared<CoreMBPTCalculator>(orbitals, bare_one_body_integrals, bare_two_body_integrals, fermi_orbitals);
    pValenceMBPTCalculator val_mbpt = std::make_shared<ValenceMBPTCalculator>(orbitals, bare_one_body_integrals, bare_two_body_integrals, fermi_orbitals);

    auto& valence = orbitals->valence;

    pOneElectronMBPT mbpt_integrals_one = std::make_shared<OneElectronMBPT>(orbitals, hf, core_mbpt, val_mbpt, identifier + ".one.int");
    pCoreValenceIntegralsMap mbpt_integrals_two = std::make_shared<CoreValenceIntegralsMap>(orbitals, core_mbpt, val_mbpt, identifier + ".two.int");
    
    // Can do these by hand, since we have a fixed calculation
    int max_pqn_in_core = 5;
    bool is_open_shell = false;

    bool use_box = !userInput.search("MBPT/--no-extra-box");
    bool include_core = !userInput.search("MBPT/--no-core");
    mbpt_integrals_one->IncludeCore(include_core, include_core && is_open_shell);
    mbpt_integrals_two->IncludeCore(include_core, include_core && is_open_shell, include_core && use_box);

    bool include_valence = userInput.search("MBPT/--use-valence");
    mbpt_integrals_one->IncludeValence(include_valence && is_open_shell);
    mbpt_integrals_two->IncludeValence(include_valence, include_valence && is_open_shell, include_valence && use_box);

    // Adjust delta
    double delta = 0.0;
    core_mbpt->SetEnergyShift(delta);
    val_mbpt->SetEnergyShift(delta);

    // Also set a floor for the valence-MBPT energy denominators
    double energy_denom_floor = 0.1;
    core_mbpt->SetEnergyFloor(energy_denom_floor);
    val_mbpt->SetEnergyFloor(energy_denom_floor);

    // Get the maximum pqn to include in two-body integrals
    std::vector<pOrbitalMap> valence_subset(4, valence);
    int max_pqn = max_pqn_in_core + 2;
    valence_subset[0].reset(new OrbitalMap(*valence));
    auto it = valence_subset[0]->begin();
    while(it != valence_subset[0]->end())
    {
        if(it->first.PQN() > max_pqn)
            it = valence_subset[0]->erase(it);
        else
            it++;
    }
    valence_subset[1] = valence_subset[0]; 

    // Now we can finally calculate the MBPT integrals. 
    // TODO: this currently needs to write the integrals to disk, which is a bit messy. Should at least delete the files when this test is done
    unsigned int one_body_size = mbpt_integrals_one->CalculateOneElectronIntegrals(valence, valence);
    unsigned int two_body_size = mbpt_integrals_two->CalculateTwoElectronIntegrals(valence_subset[0], valence_subset[1], valence_subset[2], valence_subset[3]);
    EXPECT_EQ(1, one_body_size);
    EXPECT_EQ(1, two_body_size);

    // Generate integrals
    pHFIntegrals hf_electron(new HFIntegrals(orbitals, hf));
    hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

    pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
    integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);

    // Do CI-only first
    ConfigGenerator config_generator(orbitals, userInput);
    pRelativisticConfigList relconfigs;
    Symmetry sym(1, Parity::even);
    double CI_energy = 0.;

    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();

    // Generate matrix and configurations
    auto configs = config_generator.GenerateConfigurations();
    relconfigs = config_generator.GenerateRelativisticConfigurations(configs, sym, angular_library);
    pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
    HamiltonianMatrix H_CI(hf_electron, twobody_electron, relconfigs);
    H_CI.GenerateMatrix();

    // Solve matrix
    LevelVector levels = H_CI.SolveMatrix(std::make_shared<HamiltonianID>(sym), 1);
    levels.Print();
    CI_energy = levels.levels[0]->GetEnergy();

    // Now incorporate the MBPT integrals and re-calculate
    hf_electron->Read(identifier + ".one.int");
    integrals->Read(identifier + ".two.int");

    double MBPT_energy = 0.;

    angular_library = std::make_shared<AngularDataLibrary>();

    // Generate matrix and configurations
    configs = config_generator.GenerateConfigurations();
    relconfigs = config_generator.GenerateRelativisticConfigurations(configs, sym, angular_library);
    twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
    HamiltonianMatrix H_MBPT(hf_electron, twobody_electron, relconfigs);
    H_MBPT.GenerateMatrix();

    // Solve matrix
    levels = H_MBPT.SolveMatrix(std::make_shared<HamiltonianID>(sym), 1);
    levels.Print();
    MBPT_energy = levels.levels[0]->GetEnergy();

    double energy_diff = MBPT_energy - CI_energy;
    EXPECT_NEAR(-0.01766, energy_diff, 5e-4); // TODO: check the tolerance we want to use
}

/* Comparison with valence-MBPT corrections in Blundell et al, 1988 - "Relativistic All Order Equations for Helium" .
   This tests valence-valence MBPT (no CI) in a well-characterised system.*/
TEST(CoreValenceIntegralsTester, HeCoulombPotential)
{
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);

    pLattice lattice(new Lattice(3000, 1.e-6, 50.));

    // HeI - we want to set up the MBPT options, even though we'll only use them for one component
    std::string user_input_string = std::string() + 
        "ID=HeI\n" + 
        "Z=2\n" + 
        "-s12\n" + 
        "[Lattice]\n" + 
        "NumPoints = 3000\n" + 
        "[HF]\n" + 
        "N=0\n" + 
        "Configuration=0\n" + 
        "[Basis]\n" + 
        "--bspline-basis\n" + 
        "ValenceBasis=1s\n" + 
        "[Basis/BSpline]\n" + 
        "K=9\n" + 
        "SplineType='Johnson'\n" + 
        "[CI]\n" + 
        "LeadingConfigurations='1s2'\n" + 
        "ElectronExcitations=2\n" + 
        "EvenParityTwoJ = '0'\n" + 
        "NumSolutions = 1\n" + 
        "[MBPT]\n" + 
        "--use-valence\n" + 
        "Basis=35spdfghi\n"; 
    
    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");
    std::string identifier = "He_gtest";

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();

    // Generate MBPT integrals
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pCoulombOperator coulomb(new CoulombOperator(lattice));
    pHartreeY hartreeY(new HartreeY(hf->GetIntegrator(), coulomb));
    // Bare integrals for MBPT
    pHFIntegrals bare_one_body_integrals = std::make_shared<HFIntegrals>(orbitals, hf);
    pSlaterIntegrals bare_two_body_integrals = std::make_shared<SlaterIntegralsFlatHash>(orbitals, hartreeY);

    // MBPT calculators
    std::string fermi_orbitals = "";
    pCoreMBPTCalculator core_mbpt = std::make_shared<CoreMBPTCalculator>(orbitals, bare_one_body_integrals, bare_two_body_integrals, fermi_orbitals);
    pValenceMBPTCalculator val_mbpt = std::make_shared<ValenceMBPTCalculator>(orbitals, bare_one_body_integrals, bare_two_body_integrals, fermi_orbitals);

    auto& valence = orbitals->valence;

    pOneElectronMBPT mbpt_integrals_one = std::make_shared<OneElectronMBPT>(orbitals, hf, core_mbpt, val_mbpt, identifier + ".one.int");
    pCoreValenceIntegralsMap mbpt_integrals_two = std::make_shared<CoreValenceIntegralsMap>(orbitals, core_mbpt, val_mbpt, identifier + ".two.int");
    
    // Can do these by hand, since we have a fixed calculation
    int max_pqn_in_core = 5;
    bool is_open_shell = false;

    bool use_box = !userInput.search("MBPT/--no-extra-box");
    bool include_core = !userInput.search("MBPT/--no-core");
    mbpt_integrals_one->IncludeCore(include_core, include_core && is_open_shell);
    mbpt_integrals_two->IncludeCore(include_core, include_core && is_open_shell, include_core && use_box);

    bool include_valence = userInput.search("MBPT/--use-valence");
    mbpt_integrals_one->IncludeValence(include_valence && is_open_shell);
    mbpt_integrals_two->IncludeValence(include_valence, include_valence && is_open_shell, include_valence && use_box);

    // Adjust delta
    double delta = 0.0;
    core_mbpt->SetEnergyShift(delta);
    val_mbpt->SetEnergyShift(delta);

    // Also set a floor for the valence-MBPT energy denominators
    double energy_denom_floor = 0.1;
    core_mbpt->SetEnergyFloor(energy_denom_floor);
    val_mbpt->SetEnergyFloor(energy_denom_floor);

    // Get the maximum pqn to include in two-body integrals
    std::vector<pOrbitalMap> valence_subset(4, valence);
    int max_pqn = max_pqn_in_core + 2;
    valence_subset[0].reset(new OrbitalMap(*valence));
    auto it = valence_subset[0]->begin();
    while(it != valence_subset[0]->end())
    {
        if(it->first.PQN() > max_pqn)
            it = valence_subset[0]->erase(it);
        else
            it++;
    }
    valence_subset[1] = valence_subset[0]; 

    // Now we can finally calculate the MBPT integrals. 
    // TODO: this currently needs to write the integrals to disk, which is a bit messy. Should at least delete the files when this test is done
    unsigned int one_body_size = mbpt_integrals_one->CalculateOneElectronIntegrals(valence, valence);
    unsigned int two_body_size = mbpt_integrals_two->CalculateTwoElectronIntegrals(valence_subset[0], valence_subset[1], valence_subset[2], valence_subset[3]);
    EXPECT_EQ(1, one_body_size);
    EXPECT_EQ(1, two_body_size);

    // Generate integrals
    pHFIntegrals hf_electron(new HFIntegrals(orbitals, hf));
    hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

    pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
    integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);

    // Do CI-only first
    ConfigGenerator config_generator(orbitals, userInput);
    pRelativisticConfigList relconfigs;
    Symmetry sym(0, Parity::even);
    double CI_energy = 0.;

    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();

    // Generate matrix and configurations
    auto configs = config_generator.GenerateConfigurations();
    relconfigs = config_generator.GenerateRelativisticConfigurations(configs, sym, angular_library);
    pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
    HamiltonianMatrix H_CI(hf_electron, twobody_electron, relconfigs);
    H_CI.GenerateMatrix(1); // We only have one config in the matrix

    // Solve matrix
    LevelVector levels = H_CI.SolveMatrix(std::make_shared<HamiltonianID>(sym), 1);
    levels.Print();
    CI_energy = levels.levels[0]->GetEnergy();

    // Now incorporate the MBPT integrals and re-calculate
    hf_electron->Read(identifier + ".one.int");
    integrals->Read(identifier + ".two.int");

    double MBPT_energy = 0.;

    angular_library = std::make_shared<AngularDataLibrary>();

    // Generate matrix and configurations
    configs = config_generator.GenerateConfigurations();
    relconfigs = config_generator.GenerateRelativisticConfigurations(configs, sym, angular_library);
    twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
    HamiltonianMatrix H_MBPT(hf_electron, twobody_electron, relconfigs);
    H_MBPT.GenerateMatrix(1);

    // Solve matrix
    levels = H_MBPT.SolveMatrix(std::make_shared<HamiltonianID>(sym), 1);
    levels.Print();
    MBPT_energy = levels.levels[0]->GetEnergy();

    double energy_diff = MBPT_energy - CI_energy;
    EXPECT_NEAR(-0.1576, energy_diff, 5e-4); // TODO: check the tolerance
}

