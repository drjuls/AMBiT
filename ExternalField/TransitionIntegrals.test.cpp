#include "gtest/gtest.h"
#include "Include.h"
#include "TransitionIntegrals.h"
#include "EJOperator.h"
#include "Configuration/HamiltonianMatrix.h"
#include "Basis/BasisGenerator.h"
#include "Configuration/ConfigGenerator.h"
#include "Configuration/GFactor.h"

TEST(TransitionIntegralTester, HeTransitions)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // He line strengths - comparison with Walter Johnson book page 229
    std::string user_input_string = std::string() +
        "NuclearRadius = 1.5\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 2\n" +
        "[HF]\n" +
        "N = 0\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 14spd\n" +
        "BSpline/Rmax = 50.0\n" +
        "[CI]\n" +
        "LeadingConfigurations = '1s2'\n" +
        "ElectronExcitations = 2\n";//'1,20spd,2,5spd'\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    // Generate integrals
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pOneElectronIntegrals hf_electron(new OneElectronIntegrals(orbitals, hf));
    hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

    pCoulombOperator coulomb(new CoulombOperator(lattice));
    pHartreeY hartreeY(new HartreeY(hf->GetIntegrator(), coulomb));
    pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
    integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);

    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
    ConfigGenerator config_generator(orbitals, userInput);
    pRelativisticConfigList relconfigs;
    Symmetry sym(0, Parity::even);
    double ground_state_energy = 0., excited_state_energy = 0.;

    // Generate matrix and configurations
    relconfigs = config_generator.GenerateRelativisticConfigurations(sym, angular_library);
    pTwoElectronCoulombOperator twobody_electron(new TwoElectronCoulombOperator<pSlaterIntegrals>(integrals));
    HamiltonianMatrix H_even(hf_electron, twobody_electron, relconfigs);
    H_even.GenerateMatrix();

    // Solve matrix
    auto hID = std::make_shared<HamiltonianID>(sym, relconfigs);
    LevelVector levels = H_even.SolveMatrix(hID, 3);
    Print(levels);
    pLevel ground_state = levels[0];
    pLevel singlet_S = levels[1];
    ground_state_energy = ground_state->GetEnergy();

    // Rinse and repeat for excited state
    sym = Symmetry(2, Parity::odd);
    relconfigs = config_generator.GenerateRelativisticConfigurations(sym, angular_library);
    hID = std::make_shared<HamiltonianID>(sym, relconfigs);

    HamiltonianMatrix H_odd(hf_electron, twobody_electron, relconfigs);
    H_odd.GenerateMatrix();
    levels = H_odd.SolveMatrix(hID, 3);
    GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);
    g_factors.CalculateGFactors(levels);
    Print(levels);

    pLevel singlet_P = levels[1];
    pLevel triplet_P = levels[0];
    excited_state_energy = singlet_P->GetEnergy();
    EXPECT_NEAR(1.5, triplet_P->GetgFactor(), 0.001);
    EXPECT_NEAR(1.0, singlet_P->GetgFactor(), 0.001);

    // Check energy 1s2 -> 1s2p 3P1 (should be within 10%)
    EXPECT_NEAR(171160, (excited_state_energy - ground_state_energy) * MathConstant::Instance()->HartreeEnergyInInvCm(), 17000);

    // Calculate transition
    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pSpinorMatrixElement E1(new EJOperator(constants, 1, integrator, TransitionGauge::Length));
    pTransitionIntegrals E1_matrix_elements(new TransitionIntegrals(orbitals, E1));
    E1_matrix_elements->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
    ManyBodyOperator<pTransitionIntegrals> E1_many_body(E1_matrix_elements);

    double matrix_element = E1_many_body.GetMatrixElement(*ground_state, *singlet_P);
    double reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(2, 0, 1, 2, 0);
    *logstream << "1s2 -> 2s2p 1P1: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    matrix_element = E1_many_body.GetMatrixElement(*singlet_P, *ground_state);
    reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(0, 2, 1, 0, -2);
    *logstream << "2s2p 1P1 -> 1s2: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    EXPECT_NEAR(0.5311, reduced_matrix_element * reduced_matrix_element, 0.05);

    matrix_element = E1_many_body.GetMatrixElement(*singlet_S, *singlet_P);
    reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(0, 2, 1, 0, -2);
    *logstream << "1s2s 1S0 -> 2s2p 1P1: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    matrix_element = E1_many_body.GetMatrixElement(*singlet_P, *singlet_S);
    reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(0, 2, 1, 0, -2);
    *logstream << "2s2p 1P1 -> 1s2s 1S0: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    EXPECT_NEAR(25.52, reduced_matrix_element * reduced_matrix_element, 2.5);

    matrix_element = E1_many_body.GetMatrixElement(*triplet_P, *ground_state);
    reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(2, 0, 1, 2, 0);
    *logstream << "1s2p 3P1 -> 1s2 1S0: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    matrix_element = E1_many_body.GetMatrixElement(*ground_state, *triplet_P);
    reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(2, 0, 1, 2, 0);
    *logstream << "1s2 1S0 -> 1s2p 3P1: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    EXPECT_NEAR(5.4e-8, reduced_matrix_element * reduced_matrix_element, 1.e-4);
}
