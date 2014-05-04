#include "Basis/BasisGenerator.h"
#include "gtest/gtest.h"
#include "Include.h"

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

        userInput = new MultirunOptions("template.input", "//", "\n", ",");
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
    pHFOperatorConst hf = core_generator->GetHFOperator();
    pOPIntegrator integrator = hf->GetOPIntegrator();

    // Check that < 4s | h | 4s > = h.Energy()
    pOrbitalConst core_4s = core->GetState(OrbitalInfo(4, -1));
    ASSERT_FALSE(core_4s == NULL);
    EXPECT_NEAR(1.0, core_4s->Norm(integrator), 1.e-8);
    EXPECT_NEAR(-1.333798, core->GetState(OrbitalInfo(3, -2))->Energy(), 1.e-6 * 1.33380);
    EXPECT_NEAR(core_4s->Energy(), hf->GetMatrixElement(*core_4s, *core_4s), 1.e-6 * fabs(core_4s->Energy()));
}

TEST_F(BasisGeneratorTester, BSplineBasis)
{
    DebugOptions.OutputHFExcited(true);
    pHFOperatorConst hf = core_generator->GetHFOperator();
    pOPIntegrator integrator = hf->GetOPIntegrator();
    pOrbitalMapConst excited = core_generator->GenerateBasis();

    // Check that < 4s | h | 4s > = h.Energy()
    pOrbitalConst valence_state = excited->GetState(OrbitalInfo(4, -1));
    ASSERT_FALSE(valence_state == NULL);
    EXPECT_NEAR(1.0, valence_state->Norm(integrator), 1.e-8);
    EXPECT_NEAR(-0.19623375, valence_state->Energy(), 1.e-3 * 0.196);

    EXPECT_NEAR(0.037042717, excited->GetState(OrbitalInfo(5, 2))->Energy(), 1.e-3 * 0.037);   // 5d
    EXPECT_NEAR(0.026788562, excited->GetState(OrbitalInfo(5, -4))->Energy(), 1.e-3 * 0.0267);  // 5f*
    EXPECT_NEAR(0.022253172, excited->GetState(OrbitalInfo(6, 1))->Energy(), 1.e-3 * 0.0222);   // 6p
}
