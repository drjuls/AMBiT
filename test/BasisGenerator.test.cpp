#include "Basis/BasisGenerator.h"
#include "gtest/gtest.h"
#include "Include.h"
#include "Specification/Specification.h"

using namespace Ambit;

class BasisGeneratorTester : public ::testing::Test
{
protected:
    static void SetUpTestCase() {
        lattice = pLattice(new Lattice(1000, 1.e-6, 50.));

        DebugOptions.LogFirstBuild(true);
        DebugOptions.LogHFIterations(true);

        ImportSpec();
        core_generator = new BasisGenerator(lattice);
        core = core_generator->GenerateHFCore();
    }

    static void TearDownTestCase() {
        delete core_generator;
        core_generator = NULL;
    }

    // Re-import the suite spec after each per-test Reset().
    void SetUp() override { ImportSpec(); }

    static void ImportSpec() {
        auto specification = GlobalSpecification::Instance();
        std::string perr = importSpecificationINIString(specification, user_input_strings);
        ASSERT_TRUE(perr.empty()) << "importSpecificationINIString:\n" << perr;
        perr = specification->validateAndNormaliseSpecification();
        ASSERT_TRUE(perr.empty()) << "validateAndNormaliseSpecification:\n" << perr;
    }

    static const std::string user_input_strings;
    static pLattice lattice;
    static pCore core;
    static BasisGenerator* core_generator;
};

const std::string BasisGeneratorTester::user_input_strings = std::string() +
    "-c\n" +
    "ID = test \n" +
    "NuclearRadius = 3.7188 \n" +
    "NuclearThickness = 2.3 \n" +
    "Z = 20 \n" +
    "[Lattice] \n" +
    "NumPoints = 1000 \n" +
    "StartPoint = 1.e-6 \n" +
    "EndPoint = 50. \n" +
    "[HF] \n" +
    "N = 18 \n" +
    "Configuration = '1s2 2s2 2p6 3s2 3p6' \n" +
    "[Basis] \n" +
    "--bspline-basis\n" +
    "ValenceBasis = 6spdf \n" +
    "FrozenCore = 3s2p \n" +
    "BSpline/Rmax = 50.0 \n" +
    "[CI] \n" +
    "LeadingConfigurations = '4s2, 4s1 4p1' \n" +
    "ExtraConfigurations   = '4p2' \n" +
    "ElectronExcitations = 2 \n" +
    "HoleExcitations = 1 \n" +
    "EvenParityTwoJ = '0' \n" +
    "OddParityTwoJ =  '0, 2' \n" +
    "NumSolutions = 3";
pLattice BasisGeneratorTester::lattice = pLattice();
pCore BasisGeneratorTester::core = pCore();
BasisGenerator* BasisGeneratorTester::core_generator = NULL;

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
    // Turn off BSplines and activate HF basis
    auto specification = GlobalSpecification::Instance();
    specification->basis_bspline = false;
    specification->basis_hf = true;

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
