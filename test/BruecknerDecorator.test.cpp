#include "MBPT/BruecknerDecorator.h"
#include "gtest/gtest.h"
#include "Include.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ODESolver.h"
#include "Universal/PhysicalConstant.h"
#include "HartreeFock/HartreeFocker.h"
#include "Basis/BasisGenerator.h"
#include "Universal/MathConstant.h"
#include "MBPT/CoreMBPTCalculator.h"
#include "Specification/Specification.h"

using namespace Ambit;

TEST(BruecknerDecoratorTester, MgIISlow)
{
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);

    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // MgII
    std::string user_input_string = std::string() +
        "ID = MgII\n" +
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

    auto specification = GlobalSpecification::Instance();
    // Can parse the user_input_string directly via parapara
    std::string perr = importSpecificationINIString(specification, user_input_string);
    ASSERT_TRUE(perr.empty()) << "importSpecificationINIString:\n" << perr;
    perr = specification->validateAndNormaliseSpecification();
    ASSERT_TRUE(perr.empty()) << "validateAndNormaliseSpecification:\n" << perr;
    // Get core and excited basis
    BasisGenerator basis_generator(lattice);
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
    pSlaterIntegrals two_body_integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
    CoreMBPTCalculator mbpt(orbitals, one_body_integrals, two_body_integrals, specification->mbpt_energy_denom_orbitals);
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
