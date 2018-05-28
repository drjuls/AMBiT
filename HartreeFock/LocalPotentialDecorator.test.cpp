#include "gtest/gtest.h"
#include "Include.h"
#include "LocalPotentialDecorator.h"
#include "Atom/MultirunOptions.h"
#include "Basis/BasisGenerator.h"

/** Comparisons in this section made with
        Sapirstein & Cheng, PRA 66, 042501 (2002).
 */

using namespace Ambit;

TEST(LocalPotentialDecoratorTester, CoreHartree)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Na
    std::string user_input_string = std::string() +
    "NuclearRadius = 3.0\n" +
    "NuclearThickness = 2.3\n" +
    "Z = 11\n" +
    "[HF]\n" +
    "--local-exchange\n" +
    "Xalpha = 0.0\n" +
    "N = 10\n" +
    "Configuration = '1s2 2s2 2p6'\n" +
    "[Basis]\n" +
    "--hf-basis\n" +
    "ValenceBasis = 3s\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();

    pOrbitalMapConst valence = orbitals->GetOrbitalMap(OrbitalClassification::valence);
    pOrbitalConst threeS = valence->GetState(OrbitalInfo(3, -1));

    EXPECT_NEAR(-0.173341, threeS->Energy(), 0.000002);
}

TEST(LocalPotentialDecoratorTester, DiracHartree)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Na
    std::string user_input_string = std::string() +
    "NuclearRadius = 3.0\n" +
    "NuclearThickness = 2.3\n" +
    "Z = 11\n" +
    "[HF]\n" +
    "--local-exchange\n" +
    "Xalpha = 0.0\n" +
    "N = 11\n" +
    "Configuration = '1s2 2s2 2p6: 3s1'\n" +
    "[Basis]\n" +
    "--hf-basis\n" +
    "ValenceBasis = 3s\n" +
    "BSpline/Rmax = 40.0\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOrbitalMapConst valence = orbitals->GetOrbitalMap(OrbitalClassification::valence);
    pOrbitalConst threeS = valence->GetState(OrbitalInfo(3, -1));

    EXPECT_NEAR(-0.167970, threeS->Energy(), 0.000002);
}

TEST(LocalPotentialDecoratorTester, KohnSham)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // K
    std::string user_input_string = std::string() +
    "NuclearRadius = 3.0\n" +
    "NuclearThickness = 2.3\n" +
    "Z = 19\n" +
    "[HF]\n" +
    "--local-exchange\n" +
    "Xalpha = 0.6666666666666667\n" +
    "N = 19\n" +
    "Configuration = '1s2 2s2 2p6 3s2 3p6: 4s1'\n" +
    "[Basis]\n" +
    "--hf-basis\n" +
    "ValenceBasis = 4s\n" +
    "BSpline/Rmax = 40.0\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOrbitalMapConst valence = orbitals->GetOrbitalMap(OrbitalClassification::valence);
    pOrbitalConst fourS = valence->GetState(OrbitalInfo(4, -1));

    EXPECT_NEAR(-0.144032, fourS->Energy(), 0.000002);
}

TEST(LocalPotentialDecoratorTester, DiracSlater)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Cs
    std::string user_input_string = std::string() +
    "NuclearRadius = 6.0\n" +
    "NuclearThickness = 2.3\n" +
    "Z = 55\n" +
    "[HF]\n" +
    "--local-exchange\n" +
    "Xalpha = 1.0\n" +
    "N = 55\n" +
    "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6: 6s1'\n" +
    "[Basis]\n" +
    "--hf-basis\n" +
    "ValenceBasis = 6s\n" +
    "BSpline/Rmax = 40.0\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOrbitalMapConst valence = orbitals->GetOrbitalMap(OrbitalClassification::valence);
    pOrbitalConst sixS = valence->GetState(OrbitalInfo(6, -1));

    EXPECT_NEAR(-0.135835, sixS->Energy(), 0.000002);
}
