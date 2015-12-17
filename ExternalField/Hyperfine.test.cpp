#include "gtest/gtest.h"
#include "Include.h"
#include "Hyperfine.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Basis/BasisGenerator.h"
#include "Atom/MultirunOptions.h"
#include "MBPT/OneElectronIntegrals.h"
#include "MBPT/SlaterIntegrals.h"

TEST(HyperfineTester, Rb)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Rb - comparison with Walter Johnson book
    std::string user_input_string = std::string() +
        "NuclearRadius = 4.8708\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 37\n" +
        "[HF]\n" +
        "N = 36\n" +
        "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6'\n" +
        "[Basis]\n" +
        "--hf-basis\n" +
        "ValenceBasis = 5s\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    HyperfineDipoleOperator HFS(lattice, integrator);

    const Orbital& s = *orbitals->valence->GetState(OrbitalInfo(5, -1));
    double g_I = 0.54121;
    double MHz = MathConstant::Instance()->AtomicFrequencyMHz();

    EXPECT_NEAR(643.9, HFS.GetMatrixElement(s, s) * g_I * MHz, 0.1);
}
