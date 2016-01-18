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
#include "MBPT/BruecknerDecorator.h"
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

    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    HyperfineDipoleOperator HFS(integrator);

    MathConstant* math = MathConstant::Instance();
    const Orbital& s = *orbitals->valence->GetState(OrbitalInfo(5, -1));
    double g_I = 0.54121;
    double MHz = math->AtomicFrequencyMHz();

    // Comparisons made with stretched states m = j
    EXPECT_NEAR(643.9, HFS.GetMatrixElement(s, s) * g_I/s.J() * MHz, 0.1);
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
        "--bspline-basis\n" +
        "ValenceBasis = 3spd\n" +
        "[MBPT]\n" +
        "Basis = 20spdf\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pSpinorOperator HFS = std::make_shared<HyperfineDipoleOperator>(integrator);

    MathConstant* math = MathConstant::Instance();
    pOrbital s = orbitals->valence->GetState(OrbitalInfo(3, -1));
    double g_I = 2.2176/1.5;
    double MHz = math->AtomicFrequencyMHz();

    EXPECT_NEAR(623.64, HFS->GetMatrixElement(*s, *s) * g_I/s->J() * MHz, 0.5);

    // p-wave
    s = orbitals->valence->GetState(OrbitalInfo(3, 1));
    EXPECT_NEAR(63.43, HFS->GetMatrixElement(*s, *s) * g_I/s->J() * MHz, 0.01);

    s = orbitals->valence->GetState(OrbitalInfo(3, -2));
    EXPECT_NEAR(12.60, HFS->GetMatrixElement(*s, *s) * g_I/s->J() * MHz, 0.01);

    // Get HF+HFS operator
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pRPAOperator rpa = std::make_shared<RPAOperator<HyperfineDipoleOperator>>(HFS, basis_generator.GetHartreeY());
    RPASolver rpa_solver(core);

    DebugOptions.LogHFIterations(true);
    rpa_solver.SolveRPACore(hf, rpa);

    s = orbitals->valence->GetState(OrbitalInfo(3, -1));
    EXPECT_NEAR(767.24, rpa->GetMatrixElement(*s, *s) * g_I/s->J() * MHz, 0.5);

    s = orbitals->valence->GetState(OrbitalInfo(3, 1));
    EXPECT_NEAR(82.33, rpa->GetMatrixElement(*s, *s) * g_I/s->J() * MHz, 0.01);

    s = orbitals->valence->GetState(OrbitalInfo(3, -2));
    EXPECT_NEAR(18.01, rpa->GetMatrixElement(*s, *s) * g_I/s->J() * MHz, 0.01);

    // Brueckner orbitals
    s = orbitals->valence->GetState(OrbitalInfo(3, -1));

    pOneElectronIntegrals one_body_integrals(new OneElectronIntegrals(orbitals, hf));
    pSlaterIntegrals two_body_integrals(new SlaterIntegralsMap(orbitals, basis_generator.GetHartreeY()));
    CoreMBPTCalculator mbpt(orbitals, one_body_integrals, two_body_integrals);
    mbpt.UpdateIntegrals();

    *logstream << std::setprecision(10);

    double hf_energy = s->Energy();
    *logstream << "HF:         " << hf_energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;

    double direct_summation = hf_energy + mbpt.GetOneElectronDiagrams(OrbitalInfo(3, -1), OrbitalInfo(3, -1));
    *logstream << "Direct sum: " << direct_summation * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;

    // Calculate sigma potential
    pBruecknerDecorator brueckner(new BruecknerDecorator(hf));
    brueckner->IncludeLower(false);

    pOrbital brueckner_target(new Orbital(*s));
    brueckner->CalculateSigma(brueckner_target->Kappa(), orbitals, basis_generator.GetHartreeY(), hf);

    double brueckner_matrix_element = brueckner->GetMatrixElement(*brueckner_target, *brueckner_target);
    *logstream << "Brueckner:  " << brueckner_matrix_element * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;

    // Iterate
    pODESolver ode_solver(new AdamsSolver(integrator));
    HartreeFocker hartree_focker(ode_solver);
    hartree_focker.ConvergeOrbitalAndExchange(brueckner_target, brueckner);
    *logstream << "Iterated:   " << brueckner_target->Energy() * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;

    // Comparison with experiment
    EXPECT_NEAR(885.81, rpa->GetMatrixElement(*brueckner_target, *brueckner_target) * g_I/s->J() * MHz, 2.0);
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

    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    double Rmag = 5.6748; // Nuclear magnetic radius (fm)
    pSpinorOperator HFS = std::make_shared<HyperfineDipoleOperator>(integrator, Rmag);

    MathConstant* math = MathConstant::Instance();
    double g_I = 2.582768/3.5;
    double kCm = math->HartreeEnergyInInvCm() * 1000.;

    pOrbital s, d;
    double stretched_state;     // Convert from reduced matrix element to stretched state m = j

    s = orbitals->valence->GetState(OrbitalInfo(6, -1));
    stretched_state = math->Electron3j(s->TwoJ(), s->TwoJ(), 1, -s->TwoJ(), s->TwoJ());
    EXPECT_NEAR(47.54, stretched_state * HFS->GetReducedMatrixElement(*s, *s) * g_I/s->J() * kCm, 0.01);

    d = orbitals->valence->GetState(OrbitalInfo(5, 2));
    EXPECT_NEAR(0.505e-10, HFS->GetReducedMatrixElement(*d, *s) * g_I, 0.001e-10);

    // Get HF+HFS operator
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pRPAOperator rpa = std::make_shared<RPAOperator<HyperfineDipoleOperator>>(HFS, basis_generator.GetHartreeY());
    RPASolver rpa_solver(core);

    DebugOptions.LogHFIterations(true);
    rpa_solver.SolveRPACore(hf, rpa);

    pRPAOrbital ext = std::make_shared<RPAOrbital>(*s);
    double deltaE = rpa_solver.CalculateRPAExcited(ext, rpa);
    EXPECT_NEAR(57.31, stretched_state * deltaE * g_I/s->J() * kCm, 0.01);
}
