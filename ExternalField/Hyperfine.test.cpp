#include "gtest/gtest.h"
#include "Include.h"
#include "Hyperfine.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ConfigurationParser.h"
#include "HartreeFock/HartreeFocker.h"
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

    MathConstant* math = MathConstant::Instance();
    const Orbital& s = *orbitals->valence->GetState(OrbitalInfo(5, -1));
    double g_I = 0.54121;
    double MHz = math->AtomicFrequencyMHz();

    // Convert from reduced matrix element to stretched state m = j
    double stretched_state = math->Electron3j(s.TwoJ(), s.TwoJ(), 1, -s.TwoJ(), s.TwoJ());
    EXPECT_NEAR(643.9, stretched_state * HFS.GetMatrixElement(s, s) * g_I/s.J() * MHz, 0.1);
}

TEST(HyperfineTester, Na)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Na
    std::string user_input_string = std::string() +
    "NuclearRadius = 4.8708\n" +
    "NuclearThickness = 2.3\n" +
    "Z = 11\n" +
    "[HF]\n" +
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
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    HyperfineDipoleOperator HFS(lattice, integrator);

    MathConstant* math = MathConstant::Instance();
    const Orbital& s = *orbitals->valence->GetState(OrbitalInfo(3, -1));
    double g_I = 2.2176/1.5;
    double MHz = math->AtomicFrequencyMHz();

    // Convert from reduced matrix element to stretched state m = j
    double stretched_state = math->Electron3j(s.TwoJ(), s.TwoJ(), 1, -s.TwoJ(), s.TwoJ());
    EXPECT_NEAR(623.64, stretched_state * HFS.GetMatrixElement(s, s) * g_I/s.J() * MHz, 0.1);

    // RPA: replace core orbitals with RPAOrbitals
    pCore rpa_core(new Core(lattice));
    for(const auto& orb: *core)
    {
        pRPAOrbital new_orb(new RPAOrbital(*orb.second));
        for(int twoj = mmax(1, new_orb->TwoJ()-2); twoj <= new_orb->TwoJ()+2; twoj+=2)
        {
            int kappa = math->convert_to_kappa(twoj, new_orb->GetParity());
            new_orb->deltapsi.insert(std::make_pair(kappa, std::make_shared<DeltaOrbital>(kappa, new_orb)));
        }
        rpa_core->AddState(new_orb);
    }
    rpa_core->SetOccupancies(core->GetOccupancies());

    // Get HF+HFS operator
    pHFOperator hf = std::make_shared<HyperfineDipoleRPADecorator>(basis_generator.GetClosedHFOperator(), basis_generator.GetHartreeY());
    HartreeFocker HF_Solver(std::make_shared<AdamsSolver>(integrator));

    HF_Solver.SolveCore(rpa_core, hf);

    EXPECT_NEAR(core->GetState(OrbitalInfo(2, -1))->Energy(), rpa_core->GetState(OrbitalInfo(2, -1))->Energy(), 1.e-6);
}
