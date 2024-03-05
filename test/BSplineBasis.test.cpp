#include "gtest/gtest.h"
#include "Include.h"
#include "HartreeFock/Core.h"
#include "Basis/BasisGenerator.h"
#include "Atom/MultirunOptions.h"
#include "Basis/BSplineBasis.h"

using namespace Ambit;

TEST(BSplineBasisTester, Rb)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    std::string user_input_string = std::string() +
    "NuclearRadius = 4.8708\n" +
    "NuclearThickness = 2.3\n" +
    "Z = 37\n" +
    "[HF]\n" +
    "N = 36\n" +
    "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6'\n" +
    "[Basis]\n" +
    "--bspline-basis\n" +
    "ValenceBasis = 10s\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();

    int num_splines = 40;
    int k = 7;
    double rmax = 50.0;

    int kappa = 2;
    pHFOperator hf = basis_generator.GetOpenHFOperator();
    BSplineBasis basis(lattice, num_splines, k, rmax);
    pOrbitalMap excited = basis.GenerateBSplines(hf, kappa, 10);
    pOrbital s = excited->GetState(OrbitalInfo(3, kappa));
    EXPECT_NEAR(-4.87967974917, s->Energy(), 4.e-6);

    // Test s-wave and Dirac sea
    kappa = -1;
    excited = basis.GenerateCompleteBasis(hf, kappa);
    s = excited->GetState(OrbitalInfo(5, kappa));
    EXPECT_NEAR(-0.13929159, s->Energy(), 1.e-6);

    s = excited->GetState(OrbitalInfo(1, -1));
    pOrbital neg = excited->GetState(OrbitalInfo(0, -1));
    double gap = 2./hf->GetPhysicalConstant()->GetAlphaSquared() + s->Energy();   // 2mc^2 - binding energy
    EXPECT_NEAR(gap, s->Energy() - neg->Energy(), 1.e-6 * gap);
}
