#include "gtest/gtest.h"
#include "Include.h"
#include "EJOperator.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Basis/BasisGenerator.h"
#include "Atom/MultirunOptions.h"
#include "MBPT/OneElectronIntegrals.h"
#include "MBPT/SlaterIntegrals.h"

TEST(EJOperatorTester, LiTransitions)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Li - comparison with Walter Johnson book page 214
    std::string user_input_string = std::string() +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 3\n" +
        "[HF]\n" +
        "N = 2\n" +
        "Configuration = '1s2'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 2sp\n" +
        "BSpline/Rmax = 50.0\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    EJOperator E1(constants, 1, integrator);

    const Orbital& s = *orbitals->valence->GetState(OrbitalInfo(2, -1));
    const Orbital& p1 = *orbitals->valence->GetState(OrbitalInfo(2, 1));
    const Orbital& p3 = *orbitals->valence->GetState(OrbitalInfo(2, -2));

    EXPECT_NEAR(fabs(E1.GetMatrixElement(s, p1)), fabs(E1.GetMatrixElement(p1, s)), 1.e-4);
    EXPECT_NEAR(3.3644, fabs(E1.GetMatrixElement(s, p1)), 0.0003);
    EXPECT_NEAR(4.7580, fabs(E1.GetMatrixElement(s, p3)), 0.0005);

    E1.SetGauge(TransitionGauge::Velocity);

    EXPECT_NEAR(fabs(E1.GetMatrixElement(s, p3)), fabs(E1.GetMatrixElement(p3, s)), 1.e-4);
    EXPECT_NEAR(3.4301, fabs(E1.GetMatrixElement(s, p1)), 0.0003);
    EXPECT_NEAR(4.8510, fabs(E1.GetMatrixElement(s, p3)), 0.0005);
}

TEST(EJOperatorTester, NaTransitions)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Na - comparison with Walter Johnson book page 214
    std::string user_input_string = std::string() +
    "NuclearRadius = 3.7188\n" +
    "NuclearThickness = 2.3\n" +
    "Z = 11\n" +
    "[HF]\n" +
    "N = 10\n" +
    "Configuration = '1s2 2s2 2p6'\n" +
    "[Basis]\n" +
    "--bspline-basis\n" +
    "ValenceBasis = 3sp\n" +
    "BSpline/Rmax = 50.0\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    EJOperator E1(constants, 1, integrator);

    const Orbital& s = *orbitals->valence->GetState(OrbitalInfo(3, -1));
    const Orbital& p1 = *orbitals->valence->GetState(OrbitalInfo(3, 1));
    const Orbital& p3 = *orbitals->valence->GetState(OrbitalInfo(3, -2));

    EXPECT_NEAR(fabs(E1.GetMatrixElement(s, p3)), fabs(E1.GetMatrixElement(p3, s)), 1.e-4);
    EXPECT_NEAR(3.6906, fabs(E1.GetMatrixElement(s, p1)), 0.0003);
    EXPECT_NEAR(5.2188, fabs(E1.GetMatrixElement(s, p3)), 0.0005);

    E1.SetGauge(TransitionGauge::Velocity);

    EXPECT_NEAR(fabs(E1.GetMatrixElement(s, p1)), fabs(E1.GetMatrixElement(p1, s)), 1.e-4);
    EXPECT_NEAR(3.6516, fabs(E1.GetMatrixElement(s, p1)), 0.0003);
    EXPECT_NEAR(5.1632, fabs(E1.GetMatrixElement(s, p3)), 0.0005);
}
