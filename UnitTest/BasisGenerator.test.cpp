#include "Basis/BasisGenerator.h"
#include "gtest/gtest.h"
#include "Include.h"

using namespace Ambit;

class BasisGeneratorTester : public ::testing::Test
{
protected:
    // Per-test-case set-up.
    // Called before the first test in this test case.
    // Can be omitted if not needed.
    static void SetUpTestCase() {
        lattice = pLattice(new Lattice(1000, 1.e-6, 50.));

        DebugOptions.LogFirstBuild(true);
        DebugOptions.LogHFIterations(true);

        std::string user_input_string = std::string() +
        "-c \n" +
        "ID = test\n" + 
        "Multirun = 'AlphaSquaredVariation'\n" +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 20\n" +
        "AlphaSquaredVariation = '-0.1, 0.0, 0.1'\n" 
        "[Lattice]\n" +
        "NumPoints = 1000\n" +
        "StartPoint = 1.e-6\n" +
        "EndPoint = 50.\n" +
        "[HF]\n" +
        "N = 18\n" +
        "Configuration = '1s2 2s2 2p6 3s2 3p6'\n" +
        "[Basis]\n" +
        "ValenceBasis = 6spdf\n" +
        "FrozenCore = 3s2p\n" +
        "BSpline/Rmax = 50.0\n" +
        "[CI]\n" +
        "LeadingConfigurations = '4s2, 4s1 4p1'\n" +
        "ExtraConfigurations   = '4p2'\n" +
        "ElectronExcitations = 2\n" +
        "HoleExcitations = 1\n" +
        "EvenParityTwoJ = '0'\n" +
        "OddParityTwoJ =  '0, 2'\n" +
        "NumSolutions = 3\n" +
        "[./Output]\n" +
        "Type = Standard\n" +
        "ShowgFactors = 1\n" +
        "ShowPercentages = 1\n" +
        "NumLeadingConfigurations = 5\n" +
        "MinimumDisplayedPercentage = 1\n" +
        "MaxDisplayedEnergy = -10\n" +
        "[Transitions]\n" +
        "[./Output]\n" +
        "ShowLifetime = 1\n" +
        "ShowProbability = 1\n";

        std::stringstream user_input_stream(user_input_string);
        userInput = new MultirunOptions(user_input_stream, "//", "\n", ",");
        userInput->SetRun(1);

        core_generator = new BasisGenerator(lattice, *userInput);
        core = core_generator->GenerateHFCore();
    }

    // Per-test-case tear-down.
    // Called after the last test in this test case.
    // Can be omitted if not needed.
    static void TearDownTestCase() {
        delete core_generator;
        core_generator = NULL;
        delete userInput;
        userInput = NULL;
    }

    // Some expensive resource shared by all tests.
    static pLattice lattice;
    static pCore core;
    static BasisGenerator* core_generator;
    static MultirunOptions* userInput;
};

pLattice BasisGeneratorTester::lattice = pLattice();
pCore BasisGeneratorTester::core = pCore();
BasisGenerator* BasisGeneratorTester::core_generator = NULL;
MultirunOptions* BasisGeneratorTester::userInput = NULL;

TEST_F(BasisGeneratorTester, StartCore)
{
    pHFOperatorConst hf = core_generator->GetOpenHFOperator();
    pIntegrator integrator = hf->GetIntegrator();

    // Check that < 3s | h | 3s > = h.Energy()
    pOrbitalConst core_3s = core->GetState(OrbitalInfo(3, -1));
    ASSERT_FALSE(core_3s == NULL);
    EXPECT_NEAR(1.0, core_3s->Norm(integrator), 1.e-8);
    EXPECT_NEAR(-1.8718448, core->GetState(OrbitalInfo(3, -2))->Energy(), 1.e-6 * 1.8718448);
    EXPECT_NEAR(core_3s->Energy(), hf->GetMatrixElement(*core_3s, *core_3s), 1.e-6 * fabs(core_3s->Energy()));
}

TEST_F(BasisGeneratorTester, BSplineBasis)
{
    DebugOptions.OutputHFExcited(true);
    pHFOperatorConst hf = core_generator->GetOpenHFOperator();
    pIntegrator integrator = hf->GetIntegrator();
    pOrbitalMapConst excited = core_generator->GenerateBasis()->excited;

    // Check that < 4s | h | 4s > = h.Energy()
    pOrbitalConst valence_state = excited->GetState(OrbitalInfo(4, -1));
    ASSERT_FALSE(valence_state == NULL);
    EXPECT_NEAR(1.0, valence_state->Norm(integrator), 1.e-8);
    EXPECT_NEAR(-0.41663136, valence_state->Energy(), 1.e-3 * 0.41663136);

    EXPECT_NEAR(-0.10135097, excited->GetState(OrbitalInfo(5, 2))->Energy(), 1.e-3 * 0.10135097);   // 5d
    EXPECT_NEAR(-0.080142532, excited->GetState(OrbitalInfo(5, -4))->Energy(), 1.e-3 * 0.080142532);  // 5f*
    EXPECT_NEAR(-0.095195645, excited->GetState(OrbitalInfo(6, 1))->Energy(), 1.e-3 * 0.095195645);   // 6p
}

TEST_F(BasisGeneratorTester, HFBasis)
{
    char* argv[2];
    argv[0] = new char[50];
    argv[0][0] = 0;
    argv[1] = new char[50];
    std::strcpy(argv[1], "Basis/--hf-basis");
    MultirunOptions lineInput(2, argv, ",");
    userInput->absorb(lineInput);

    DebugOptions.OutputHFExcited(true);
    pHFOperatorConst hf = core_generator->GetOpenHFOperator();
    pIntegrator integrator = hf->GetIntegrator();
    pOrbitalMapConst excited = core_generator->GenerateBasis()->excited;
    
    // Check that < 4s | h | 4s > = h.Energy()
    pOrbitalConst valence_state = excited->GetState(OrbitalInfo(4, -1));
    ASSERT_FALSE(valence_state == NULL);
    EXPECT_NEAR(1.0, valence_state->Norm(integrator), 1.e-8);
    EXPECT_NEAR(-0.41663136, valence_state->Energy(), 1.e-3 * 0.41663136);
    
    EXPECT_NEAR(-0.10135097, excited->GetState(OrbitalInfo(5, 2))->Energy(), 1.e-3 * 0.10135097);   // 5d
    EXPECT_NEAR(-0.080142532, excited->GetState(OrbitalInfo(5, -4))->Energy(), 1.e-3 * 0.080142532);  // 5f*
    EXPECT_NEAR(-0.095195645, excited->GetState(OrbitalInfo(6, 1))->Energy(), 1.e-3 * 0.095195645);   // 6p
}
