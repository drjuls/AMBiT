#include "gtest/gtest.h"
#include "Include.h"
#include "RadiativePotential.h"
#include "Atom/MultirunOptions.h"
#include "Basis/BasisGenerator.h"

TEST(RadiativePotentialTester, Uehling)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Na - comparison with Ginges & Berengut
    std::string user_input_string = std::string() +
        "NuclearRadius = 2.9374\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 11\n" +
        "[HF]\n" +
        "N = 10\n" +
        "Configuration = '1s2 2s2 2p6'\n" +
        "[Basis]\n" +
        "--hf-basis\n" +
        "ValenceBasis = 3sp\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pHFOperator ueh = pHFOperator(new UehlingDecorator(hf, basis_generator.GetNuclearRMSRadius()));

    pOrbitalConst s = orbitals->valence->GetState(OrbitalInfo(3, -1));

    EXPECT_NEAR(-5.559e-7, ueh->GetMatrixElement(*s, *s) - hf->GetMatrixElement(*s, *s), 1.e-10);
}
