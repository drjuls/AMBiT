#include "BruecknerDecorator.h"
#include "gtest/gtest.h"
#include "Include.h"
#include "Atom/MultirunOptions.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ODESolver.h"
#include "Universal/PhysicalConstant.h"
#include "HartreeFock/HartreeFocker.h"
#include "Basis/BasisGenerator.h"
#include "Universal/MathConstant.h"

TEST(BruecknerDecoratorTester, MgIISlow)
{
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);

    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // MgII
    std::string user_input_string = std::string() +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 12\n" +
        "[HF]\n" +
        "N = 10\n" +
        "Configuration = '1s2 2s2 2p6'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 6spdf\n" +
        "BSpline/Rmax = 45.0\n" +
        "[MBPT]\n" +
        "Basis = 10spdf\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    lattice->resize(core->LargestOrbitalSize());
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();

    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pHartreeY hartreeY = basis_generator.GetHartreeY();

    DebugOptions.OutputHFExcited(false);
    DebugOptions.LogHFIterations(false);

    // Target orbital
    OrbitalInfo target(3, -1);

    // Calculate matrix element using core-valence MBPT
    pHFIntegrals one_body_integrals(new HFIntegrals(orbitals, hf));
    pSlaterIntegrals two_body_integrals(new SlaterIntegralsMap(orbitals, hartreeY));
    CoreMBPTCalculator mbpt(orbitals, one_body_integrals, two_body_integrals);
    mbpt.UpdateIntegrals();

    pOrbitalConst bare_target = orbitals->excited->GetState(target);
    double hf_energy = bare_target->Energy();
    *logstream << "HF:         " << hf_energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;

    double direct_summation = hf_energy + mbpt.GetOneElectronDiagrams(target, target);
    *logstream << "Direct sum: " << direct_summation * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;

    // Calculate sigma potential
    pBruecknerDecorator brueckner(new BruecknerDecorator(hf));
    brueckner->IncludeLower(false);

    pOrbital brueckner_target(new Orbital(*bare_target));
    brueckner->CalculateSigma(brueckner_target->Kappa(), orbitals, hartreeY, hf);

    double brueckner_matrix_element = brueckner->GetMatrixElement(*brueckner_target, *brueckner_target);
    *logstream << "Brueckner:  " << brueckner_matrix_element * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;

    EXPECT_NEAR(direct_summation - hf_energy, brueckner_matrix_element - hf_energy, 0.01 * fabs(direct_summation - hf_energy));

    // Iterate: usually around 10% difference
    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    HartreeFocker hartree_focker(ode_solver);
    hartree_focker.CalculateExcitedState(brueckner_target, brueckner);
    *logstream << "Iterated:   " << brueckner_target->Energy() * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;

    EXPECT_NEAR(brueckner_matrix_element - hf_energy, brueckner_target->Energy() - hf_energy, 0.2 * fabs(brueckner_matrix_element - hf_energy));
}
