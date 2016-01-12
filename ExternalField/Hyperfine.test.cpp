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
#include "RPASolver.h"

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
    HyperfineDipoleOperator HFS(integrator);

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
        "NuclearRadius = 2.8853\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 11\n" +
        "[HF]\n" +
        "N = 10\n" +
        "Configuration = '1s2 2s2 2p6'\n" +
        "[Basis]\n" +
        "--hf-basis\n" +
        "ValenceBasis = 3spd\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    HyperfineDipoleOperator HFS(integrator);

    MathConstant* math = MathConstant::Instance();
    pOrbital s = orbitals->valence->GetState(OrbitalInfo(3, -1));
    double g_I = 2.2176/1.5;
    double MHz = math->AtomicFrequencyMHz();

    // Convert from reduced matrix element to stretched state m = j
    double stretched_state = math->Electron3j(s->TwoJ(), s->TwoJ(), 1, -s->TwoJ(), s->TwoJ());
    EXPECT_NEAR(623.64, stretched_state * HFS.GetMatrixElement(*s, *s) * g_I/s->J() * MHz, 0.5);

    // p-wave
    s = orbitals->valence->GetState(OrbitalInfo(3, 1));
    stretched_state = math->Electron3j(s->TwoJ(), s->TwoJ(), 1, -s->TwoJ(), s->TwoJ());
    EXPECT_NEAR(63.43, stretched_state * HFS.GetMatrixElement(*s, *s) * g_I/s->J() * MHz, 0.01);

    s = orbitals->valence->GetState(OrbitalInfo(3, -2));
    stretched_state = math->Electron3j(s->TwoJ(), s->TwoJ(), 1, -s->TwoJ(), s->TwoJ());
    EXPECT_NEAR(12.60, stretched_state * HFS.GetMatrixElement(*s, *s) * g_I/s->J() * MHz, 0.01);

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
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pHyperfineRPAOperator rpa = std::make_shared<HyperfineRPAOperator>(rpa_core, basis_generator.GetHartreeY());
    RPASolver rpa_solver;

    DebugOptions.LogHFIterations(true);
    rpa_solver.SolveRPACore(rpa_core, hf, rpa, true);

    s = orbitals->valence->GetState(OrbitalInfo(3, -1));
    stretched_state = math->Electron3j(s->TwoJ(), s->TwoJ(), 1, -s->TwoJ(), s->TwoJ());
    pRPAOrbital ext = std::make_shared<RPAOrbital>(*s);
    double deltaE = rpa_solver.CalculateRPAExcited(ext, rpa);
    EXPECT_NEAR(767.24, stretched_state * deltaE * g_I/s->J() * MHz, 0.5);

    s = orbitals->valence->GetState(OrbitalInfo(3, 1));
    stretched_state = math->Electron3j(s->TwoJ(), s->TwoJ(), 1, -s->TwoJ(), s->TwoJ());
    ext = std::make_shared<RPAOrbital>(*s);
    deltaE = rpa_solver.CalculateRPAExcited(ext, rpa);
    EXPECT_NEAR(82.33, stretched_state * deltaE * g_I/s->J() * MHz, 0.01);

    s = orbitals->valence->GetState(OrbitalInfo(3, -2));
    stretched_state = math->Electron3j(s->TwoJ(), s->TwoJ(), 1, -s->TwoJ(), s->TwoJ());
    ext = std::make_shared<RPAOrbital>(*s);
    deltaE = rpa_solver.CalculateRPAExcited(ext, rpa);
    EXPECT_NEAR(18.01, stretched_state * deltaE * g_I/s->J() * MHz, 0.01);
}

TEST(HyperfineTester, CsRPA)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Cs - Comparison with Dzuba code
    // [private communication; similar to Dzuba, Flambaum, Sushkov, JPB 17, 1953 (1984)]
    std::string user_input_string = std::string() +
    "NuclearRadius = 5.6748\n" +
    "NuclearThickness = 2.3\n" +
    "Z = 55\n" +
    "[HF]\n" +
    "N = 54\n" +
    "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6'\n" +
    "[Basis]\n" +
    "--hf-basis\n" +
    "ValenceBasis = 6spd\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    double Rmag = 5.6748; // Nuclear magnetic radius (fm)
    HyperfineDipoleOperator HFS(integrator, Rmag);

    MathConstant* math = MathConstant::Instance();
    double g_I = 2.582768/3.5;
    double kCm = math->HartreeEnergyInInvCm() * 1000.;

    pOrbital s, d;
    double stretched_state;     // Convert from reduced matrix element to stretched state m = j

    s = orbitals->valence->GetState(OrbitalInfo(6, -1));
    stretched_state = math->Electron3j(s->TwoJ(), s->TwoJ(), 1, -s->TwoJ(), s->TwoJ());
    EXPECT_NEAR(47.54, stretched_state * HFS.GetMatrixElement(*s, *s) * g_I/s->J() * kCm, 0.01);

    d = orbitals->valence->GetState(OrbitalInfo(5, 2));
    EXPECT_NEAR(0.505e-10, HFS.GetMatrixElement(*d, *s) * g_I, 0.001e-10);

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
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pHyperfineRPAOperator rpa = std::make_shared<HyperfineRPAOperator>(rpa_core, basis_generator.GetHartreeY(), 5.6748);
    RPASolver rpa_solver;

    DebugOptions.LogHFIterations(true);
    rpa_solver.SolveRPACore(rpa_core, hf, rpa, true);

    pRPAOrbital ext = std::make_shared<RPAOrbital>(*s);
    double deltaE = rpa_solver.CalculateRPAExcited(ext, rpa);
    EXPECT_NEAR(57.31, stretched_state * deltaE * g_I/s->J() * kCm, 0.01);
}
