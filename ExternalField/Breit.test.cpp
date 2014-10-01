#include "gtest/gtest.h"
#include "Include.h"
#include "BreitHFDecorator.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Basis/BasisGenerator.h"
#include "Atom/MultirunOptions.h"
#include "MBPT/OneElectronIntegrals.h"
#include "MBPT/SlaterIntegrals.h"
#include "HartreeFock/HartreeFocker.h"

TEST(BreitTester, LiLikeNe)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 20.));

    // Li-like Ne - comparison with Walter Johnson book page 208
    std::string user_input_string = std::string() +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 10\n" +
        "[HF]\n" +
        "N = 2\n" +
        "Configuration = '1s2'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 2sp\n" +
        "BSpline/Rmax = 20.0\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));

    const Orbital& s_core = *orbitals->core->GetState(OrbitalInfo(1, -1));
    const Orbital& s = *orbitals->valence->GetState(OrbitalInfo(2, -1));
    const Orbital& p1 = *orbitals->valence->GetState(OrbitalInfo(2, 1));
    const Orbital& p3 = *orbitals->valence->GetState(OrbitalInfo(2, -2));

    pHartreeY breit(new BreitZero(pHartreeY(new HartreeYBase()), integrator, coulomb));
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pHFOperator breit_hf(new BreitHFDecorator(hf, breit));

    // Check values against Johnson
    double matrix_element = breit_hf->GetMatrixElement(s, s);
    matrix_element -= hf->GetMatrixElement(s, s);
    EXPECT_NEAR(0.00090, matrix_element, 1.e-5);

    // Check that ApplyTo and GetMatrixElement are equivalent in BreitZero for N term
    breit->SetParameters(1, s_core, s);
    double value = integrator->GetInnerProduct(s, breit->ApplyTo(s_core, s.Kappa()))/2.0;
    EXPECT_NEAR(matrix_element, value, 1.e-7);

    // Check p-wave against Johnson
    matrix_element = breit_hf->GetMatrixElement(p1, p1);
    matrix_element -= hf->GetMatrixElement(p1, p1);
    EXPECT_NEAR(0.00160, matrix_element, 1.e-4);

    // Check that ApplyTo and GetMatrixElement are equivalent in BreitZero for M and O terms
    breit->SetParameters(1, s_core, p1);
    matrix_element = breit->GetMatrixElement(p1, s_core)/2.0;
    value = integrator->GetInnerProduct(p1, breit->ApplyTo(s_core, p1.Kappa()))/2.0;
    EXPECT_NEAR(matrix_element, value, 1.e-7);

    // Check p_3/2-wave against Johnson
    matrix_element = breit_hf->GetMatrixElement(p3, p3);
    matrix_element -= hf->GetMatrixElement(p3, p3);
    EXPECT_NEAR(0.00074, matrix_element, 1.e-5);

    // Check that ApplyTo and GetMatrixElement are equivalent in BreitHFDecorator
    matrix_element = breit_hf->GetMatrixElement(p3, p3);
    value = integrator->GetInnerProduct(p3, breit_hf->ApplyTo(p3));
    EXPECT_NEAR(matrix_element, value, 1.e-7);
}

TEST(BreitTester, Hg)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);

    //Hg - comparison with Lindroth et al. J. Phys. B 22, 2447 (1989)
    std::string user_input_string = std::string() +
        "NuclearRadius = 6.0\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 80\n" +
        "[HF]\n" +
        "N = 80\n" +
        "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2'\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));

    const Orbital& one_s = *orbitals->core->GetState(OrbitalInfo(1, -1));
    const Orbital& four_d5 = *orbitals->core->GetState(OrbitalInfo(4, -3));
    const Orbital& four_f7 = *orbitals->core->GetState(OrbitalInfo(4, 3));
    const Orbital& six_s = *orbitals->core->GetState(OrbitalInfo(6, -1));

    pHartreeY breit(new BreitZero(pHartreeY(new HartreeYBase()), integrator, coulomb));
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pHFOperator breit_hf(new BreitHFDecorator(hf, breit));

    double matrix_element = breit_hf->GetMatrixElement(one_s, one_s);
    matrix_element -= hf->GetMatrixElement(one_s, one_s);
    EXPECT_NEAR(11.61, matrix_element, 0.1);

    matrix_element = breit_hf->GetMatrixElement(four_d5, four_d5);
    matrix_element -= hf->GetMatrixElement(four_d5, four_d5);
    EXPECT_NEAR(0.04968, matrix_element, 0.01);

    matrix_element = breit_hf->GetMatrixElement(four_f7, four_f7);
    matrix_element -= hf->GetMatrixElement(four_f7, four_f7);
    EXPECT_NEAR(0.02200, matrix_element, 0.01);

    matrix_element = breit_hf->GetMatrixElement(six_s, six_s);
    matrix_element -= hf->GetMatrixElement(six_s, six_s);
    EXPECT_NEAR(0.001547, matrix_element, 0.001);
}

TEST(BreitTester, IrSixteenPlus)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 20.));
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);

    //Ir16+ - comparison with Borschevsky
    std::string user_input_string = std::string() +
        "NuclearRadius = 6.0\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 77\n" +
        "[HF]\n" +
        "N = 60\n" +
        "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 5sp\n" +
        "BSpline/Rmax = 20.0\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));

    const Orbital s_hf = *orbitals->valence->GetState(OrbitalInfo(5, -1));
    const Orbital p1_hf = *orbitals->valence->GetState(OrbitalInfo(5, 1));
    const Orbital p3_hf = *orbitals->valence->GetState(OrbitalInfo(5, -2));

    pHartreeY breit(new BreitZero(pHartreeY(new HartreeYBase()), integrator, coulomb));
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pHFOperator breit_hf(new BreitHFDecorator(hf, breit));

    HartreeFocker hf_solver(ode_solver);
    hf_solver.SolveCore(core, breit_hf);

    const Orbital& s = *orbitals->valence->GetState(OrbitalInfo(5, -1));
    const Orbital& p3 = *orbitals->valence->GetState(OrbitalInfo(5, -2));

    double energy = hf->GetMatrixElement(p3, p3) - hf->GetMatrixElement(s, s);
    *logstream << "s->p3: Energy = " << energy << ", breit (%) = " << ((breit_hf->GetMatrixElement(p3, p3) - breit_hf->GetMatrixElement(s, s))/energy -1.) * 100. << std::endl;
}
