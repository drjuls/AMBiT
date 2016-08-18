#include "gtest/gtest.h"
#include "Include.h"
#include "EJOperator.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Basis/BasisGenerator.h"
#include "Atom/MultirunOptions.h"
#include "MBPT/OneElectronIntegrals.h"
#include "MBPT/SlaterIntegrals.h"
#include "Configuration/HamiltonianMatrix.h"
#include "Configuration/ConfigGenerator.h"
#include "Configuration/GFactor.h"
#include "RPAOperator.h"

TEST(EJOperatorTester, LiTransitions)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Li - comparison with Walter Johnson book page 239 (Table 8.1)
    std::string user_input_string = std::string() +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 3\n" +
        "[HF]\n" +
        "N = 2\n" +
        "Configuration = '1s2'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 2sp\n" +
        "BSpline/Rmax = 50.0\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pSpinorOperator pE1 = std::make_shared<EJOperator>(1, integrator);
    EJOperator& E1 = *std::static_pointer_cast<EJOperator>(pE1);

    const Orbital& s = *orbitals->valence->GetState(OrbitalInfo(2, -1));
    const Orbital& p1 = *orbitals->valence->GetState(OrbitalInfo(2, 1));
    const Orbital& p3 = *orbitals->valence->GetState(OrbitalInfo(2, -2));

    E1.SetFrequency(s.Energy() - p1.Energy());
    EXPECT_NEAR(fabs(E1.GetReducedMatrixElement(s, p1)), fabs(E1.GetReducedMatrixElement(p1, s)), 1.e-4);
    EXPECT_NEAR(3.3644, fabs(E1.GetReducedMatrixElement(s, p1)), 0.0001);

    E1.SetFrequency(p3.Energy() - s.Energy());
    EXPECT_NEAR(4.7580, fabs(E1.GetReducedMatrixElement(s, p3)), 0.0001);

    E1.SetGauge(TransitionGauge::Velocity);

    EXPECT_NEAR(fabs(E1.GetReducedMatrixElement(s, p3)), fabs(E1.GetReducedMatrixElement(p3, s)), 1.e-4);
    EXPECT_NEAR(4.8510, fabs(E1.GetReducedMatrixElement(s, p3)), 0.0001);
    E1.SetFrequency(s.Energy() - p1.Energy());
    EXPECT_NEAR(3.4301, fabs(E1.GetReducedMatrixElement(s, p1)), 0.0001);

    // RPA - comparison with Walter Johnson book page 240 (Table 8.2)
    DebugOptions.LogHFIterations(true);
    E1.SetGauge(TransitionGauge::Length);
    double scale = 0.001;
    pRPASolver rpa_solver = std::make_shared<RPASolver>(lattice);
    pRPAOperator rpa = std::make_shared<RPAOperator>(pE1, basis_generator.GetClosedHFOperator(), basis_generator.GetHartreeY(), rpa_solver);
    rpa->SetScale(scale);
    rpa->SetFrequency(p1.Energy() - s.Energy());
    double s_p1_length = fabs(rpa->GetReducedMatrixElement(s, p1))/scale;
    EXPECT_NEAR(3.3505, s_p1_length, 0.0001);
    rpa->SetFrequency(p3.Energy() - s.Energy());
    double s_p3_length = fabs(rpa->GetReducedMatrixElement(s, p3))/scale;
    EXPECT_NEAR(4.7383, s_p3_length, 0.0002);

//    E1.SetGauge(TransitionGauge::Velocity);
//    rpa->ClearRPACore();
//    rpa->SetFrequency(p1.Energy() - s.Energy());
//    EXPECT_NEAR(s_p1_length, fabs(rpa->GetReducedMatrixElement(s, p1))/scale, 1.e-6);
//    rpa->SetFrequency(p3.Energy() - s.Energy());
//    EXPECT_NEAR(s_p3_length, fabs(rpa->GetReducedMatrixElement(s, p3))/scale, 1.e-6);
}

TEST(EJOperatorTester, NaTransitions)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Na - comparison with Walter Johnson book page 239 (Table 8.1)
    std::string user_input_string = std::string() +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 11\n" +
        "[HF]\n" +
        "N = 10\n" +
        "Configuration = '1s2 2s2 2p6'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 3sp\n" +
        "BSpline/Rmax = 50.0\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();

    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pSpinorOperator pE1 = std::make_shared<EJOperator>(1, integrator);
    EJOperator& E1 = *std::static_pointer_cast<EJOperator>(pE1);

    const Orbital& s = *orbitals->valence->GetState(OrbitalInfo(3, -1));
    const Orbital& p1 = *orbitals->valence->GetState(OrbitalInfo(3, 1));
    const Orbital& p3 = *orbitals->valence->GetState(OrbitalInfo(3, -2));

    E1.SetFrequency(s.Energy() - p3.Energy());
    EXPECT_NEAR(fabs(E1.GetReducedMatrixElement(s, p3)), fabs(E1.GetReducedMatrixElement(p3, s)), 1.e-4);
    EXPECT_NEAR(5.2188, fabs(E1.GetReducedMatrixElement(s, p3)), 0.0005);
    E1.SetFrequency(s.Energy() - p1.Energy());
    EXPECT_NEAR(3.6906, fabs(E1.GetReducedMatrixElement(s, p1)), 0.0003);

    E1.SetGauge(TransitionGauge::Velocity);

    EXPECT_NEAR(fabs(E1.GetReducedMatrixElement(s, p1)), fabs(E1.GetReducedMatrixElement(p1, s)), 1.e-4);
    EXPECT_NEAR(3.6516, fabs(E1.GetReducedMatrixElement(s, p1)), 0.0003);
    E1.SetFrequency(p3.Energy() - s.Energy());
    EXPECT_NEAR(5.1632, fabs(E1.GetReducedMatrixElement(s, p3)), 0.0005);

    // RPA - comparison with Walter Johnson book page 240 (Table 8.2)
    //  and comparison of length and velocity gauges.
    E1.SetGauge(TransitionGauge::Length);
    pRPASolver rpa_solver = std::make_shared<RPASolver>(lattice);
    pRPAOperator rpa = std::make_shared<RPAOperator>(pE1, basis_generator.GetClosedHFOperator(), basis_generator.GetHartreeY(), rpa_solver);
    rpa->SetFrequency(p1.Energy() - s.Energy());
    double s_p1_length = fabs(rpa->GetReducedMatrixElement(s, p1));
    EXPECT_NEAR(3.6474, s_p1_length, 0.0003);
    rpa->SetFrequency(p3.Energy() - s.Energy());
    double s_p3_length = fabs(rpa->GetReducedMatrixElement(s, p3));
    EXPECT_NEAR(5.1578, s_p3_length, 0.0004);

//    DebugOptions.LogHFIterations(true);
//    E1.SetGauge(TransitionGauge::Velocity);
//    rpa->SetFrequency(p1.Energy() - s.Energy());
//    EXPECT_NEAR(s_p1_length, fabs(rpa->GetReducedMatrixElement(s, p1)), 1.e-6);
//    rpa->SetFrequency(p3.Energy() - s.Energy());
//    EXPECT_NEAR(s_p3_length, fabs(rpa->GetReducedMatrixElement(s, p3)), 1.e-6);
}

TEST(MJOperatorTester, LiTransitions)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Li fine-structure - comparison with Walter Johnson book page 163
    std::string user_input_string = std::string() +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 3\n" +
        "[HF]\n" +
        "N = 2\n" +
        "Configuration = '1s2'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 2sp\n" +
        "BSpline/Rmax = 50.0\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();

    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    MJOperator M1(1, integrator);

    const Orbital& p1 = *orbitals->valence->GetState(OrbitalInfo(2, 1));
    const Orbital& p3 = *orbitals->valence->GetState(OrbitalInfo(2, -2));

    M1.SetFrequency(p1.Energy() - p3.Energy());
    EXPECT_NEAR(fabs(M1.GetReducedMatrixElement(p3, p1)), fabs(M1.GetReducedMatrixElement(p1, p3)), 1.e-4);
    EXPECT_NEAR(2./sqrt(3.), fabs(M1.GetReducedMatrixElement(p3, p1)), 0.0003);
}

TEST(EJOperatorTester, HeTransitions)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // He line strengths - comparison with Walter Johnson book page 229
    std::string user_input_string = std::string() +
    "NuclearRadius = 1.5\n" +
    "NuclearThickness = 2.3\n" +
    "Z = 2\n" +
    "[HF]\n" +
    "N = 0\n" +
    "[Basis]\n" +
    "--bspline-basis\n" +
    "ValenceBasis = 14spd\n" +
    "BSpline/Rmax = 50.0\n" +
    "[CI]\n" +
    "LeadingConfigurations = '1s2'\n" +
    "ElectronExcitations = 2\n";//'1,20spd,2,5spd'\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    // Generate integrals
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pHFIntegrals hf_electron(new HFIntegrals(orbitals, hf));
    hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

    pCoulombOperator coulomb(new CoulombOperator(lattice));
    pHartreeY hartreeY(new HartreeY(hf->GetIntegrator(), coulomb));
    pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
    integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);

    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
    ConfigGenerator config_generator(orbitals, userInput);
    pRelativisticConfigList relconfigs;
    Symmetry sym(0, Parity::even);
    double ground_state_energy = 0., excited_state_energy = 0.;

    // Generate matrix and configurations
    relconfigs = config_generator.GenerateRelativisticConfigurations(sym, angular_library);
    pTwoElectronCoulombOperator twobody_electron(new TwoElectronCoulombOperator<pSlaterIntegrals>(integrals));
    HamiltonianMatrix H_even(hf_electron, twobody_electron, relconfigs);
    H_even.GenerateMatrix();

    // Solve matrix
    auto hID = std::make_shared<HamiltonianID>(sym, relconfigs);
    LevelVector levels = H_even.SolveMatrix(hID, 3);
    Print(levels);
    pLevel ground_state = levels[0];
    pLevel singlet_S = levels[1];
    ground_state_energy = ground_state->GetEnergy();

    // Rinse and repeat for excited state
    sym = Symmetry(2, Parity::odd);
    relconfigs = config_generator.GenerateRelativisticConfigurations(sym, angular_library);
    hID = std::make_shared<HamiltonianID>(sym, relconfigs);

    HamiltonianMatrix H_odd(hf_electron, twobody_electron, relconfigs);
    H_odd.GenerateMatrix();
    levels = H_odd.SolveMatrix(hID, 3);
    GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);
    g_factors.CalculateGFactors(levels);
    Print(levels);

    pLevel singlet_P = levels[1];
    pLevel triplet_P = levels[0];
    excited_state_energy = singlet_P->GetEnergy();
    EXPECT_NEAR(1.5, triplet_P->GetgFactor(), 0.001);
    EXPECT_NEAR(1.0, singlet_P->GetgFactor(), 0.001);

    // Check energy 1s2 -> 1s2p 3P1 (should be within 10%)
    EXPECT_NEAR(171160, (excited_state_energy - ground_state_energy) * MathConstant::Instance()->HartreeEnergyInInvCm(), 17000);

    // Calculate transition
    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pTimeDependentSpinorOperator E1(new EJOperator(1, integrator, TransitionGauge::Length));
    pTransitionIntegrals E1_matrix_elements(new TransitionIntegrals(orbitals, E1));
    ManyBodyOperator<pTransitionIntegrals> E1_many_body(E1_matrix_elements);

    E1->SetFrequency(ground_state_energy - excited_state_energy);
    E1_matrix_elements->clear();
    E1_matrix_elements->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
    double matrix_element = E1_many_body.GetMatrixElement(*ground_state, *singlet_P);
    double reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(2, 0, 1, 2, 0);
    *logstream << "1s2 -> 1s2p 1P1: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    matrix_element = E1_many_body.GetMatrixElement(*singlet_P, *ground_state);
    reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(0, 2, 1, 0, -2);
    *logstream << "1s2p 1P1 -> 1s2: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    EXPECT_NEAR(0.5311, reduced_matrix_element * reduced_matrix_element, 0.05);

    E1->SetFrequency(singlet_S->GetEnergy() - excited_state_energy);
    E1_matrix_elements->clear();
    E1_matrix_elements->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
    matrix_element = E1_many_body.GetMatrixElement(*singlet_S, *singlet_P);
    reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(0, 2, 1, 0, -2);
    *logstream << "1s2s 1S0 -> 1s2p 1P1: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    matrix_element = E1_many_body.GetMatrixElement(*singlet_P, *singlet_S);
    reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(0, 2, 1, 0, -2);
    *logstream << "1s2p 1P1 -> 1s2s 1S0: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    EXPECT_NEAR(25.52, reduced_matrix_element * reduced_matrix_element, 2.5);

    E1->SetFrequency(triplet_P->GetEnergy() - ground_state_energy);
    E1_matrix_elements->clear();
    E1_matrix_elements->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
    matrix_element = E1_many_body.GetMatrixElement(*triplet_P, *ground_state);
    reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(2, 0, 1, 2, 0);
    *logstream << "1s2p 3P1 -> 1s2 1S0: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    matrix_element = E1_many_body.GetMatrixElement(*ground_state, *triplet_P);
    reduced_matrix_element = matrix_element/MathConstant::Instance()->Electron3j(2, 0, 1, 2, 0);
    *logstream << "1s2 1S0 -> 1s2p 3P1: S_E1 = " << reduced_matrix_element * reduced_matrix_element << std::endl;
    EXPECT_NEAR(5.4e-8, reduced_matrix_element * reduced_matrix_element, 1.e-4);
}

TEST(MJOperatorTester, HolesVsElectrons)
{
    // Test holes using same CuIII CI from HamiltonianMatrixTester::HolesVsElectrons,
    // assuming that one passes
    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(false);

    pLattice lattice(new Lattice(1000, 1.e-6, 50.));
    MathConstant* math = MathConstant::Instance();
    Symmetry sym(5, Parity::even);

    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
    double m1_electron, m1_hole;

    // CuIII - old electron way
    {   // Block for reusing names
        std::string user_input_string = std::string() +
            "NuclearRadius = 3.7188\n" +
            "NuclearThickness = 2.3\n" +
            "Z = 29\n" +
            "[HF]\n" +
            "N = 28\n" +
            "Configuration = '1s2 2s2 2p6 3s2 3p6: 3d10'\n" +
            "[Basis]\n" +
            "--bspline-basis\n" +
            "ValenceBasis = 5spd\n" +
            "BSpline/Rmax = 50.0\n" +
            "[CI]\n" +
            "LeadingConfigurations = '3d9'\n" +
            "ElectronExcitations = 1\n";

        std::stringstream user_input_stream(user_input_string);
        MultirunOptions userInput(user_input_stream, "//", "\n", ",");

        // Get core and excited basis
        BasisGenerator basis_generator(lattice, userInput);
        basis_generator.GenerateHFCore();
        pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
        ConfigGenerator gen(orbitals, userInput);

        // Generate integrals
        pHFOperator hf = basis_generator.GetClosedHFOperator();
        pHFIntegrals hf_electron(new HFIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

        pCoulombOperator coulomb(new CoulombOperator(lattice));
        pHartreeY hartreeY(new HartreeY(hf->GetIntegrator(), coulomb));
        pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron(new TwoElectronCoulombOperator<pSlaterIntegrals>(integrals));

        pHamiltonianID key = std::make_shared<HamiltonianID>(sym);
        pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(sym, angular_library);
        HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
        H.GenerateMatrix();
        LevelVector levels = H.SolveMatrix(key, 3);

        // Calculate transition
        pTimeDependentSpinorOperator M1 = std::make_shared<MJOperator>(1, hf->GetIntegrator());
        M1->SetFrequency(levels[2]->GetEnergy() - levels[1]->GetEnergy());
        pTransitionIntegrals M1_matrix_elements(new TransitionIntegrals(orbitals, M1));
        M1_matrix_elements->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
        ManyBodyOperator<pTransitionIntegrals> M1_many_body(M1_matrix_elements);

        m1_electron = M1_many_body.GetMatrixElement(*levels[2], *levels[1]);

        // Convert to strength = (reduced matrix element)^2
        m1_electron = m1_electron/math->Electron3j(sym.GetTwoJ(), sym.GetTwoJ(), 1, sym.GetTwoJ(), -sym.GetTwoJ());
        m1_electron = m1_electron * m1_electron;
    }

    ASSERT_NE(m1_electron, 0.0);

    // CuIII - using holes now
    {   // Block for reusing names
        std::string user_input_string = std::string() +
            "NuclearRadius = 3.7188\n" +
            "NuclearThickness = 2.3\n" +
            "Z = 29\n" +
            "[HF]\n" +
            "N = 28\n" +
            "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10'\n" +
            "[Basis]\n" +
            "--bspline-basis\n" +
            "ValenceBasis = 5spd\n" +
            "FrozenCore = 3sp\n" +
            "BSpline/Rmax = 50.0\n" +
            "[CI]\n" +
            "LeadingConfigurations = '3d-1'\n" +
            "ElectronExcitations = 1\n" +
            "HoleExcitations = 1\n";

        std::stringstream user_input_stream(user_input_string);
        MultirunOptions userInput(user_input_stream, "//", "\n", ",");

        // Get core and excited basis
        BasisGenerator basis_generator(lattice, userInput);
        basis_generator.GenerateHFCore();
        pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
        ConfigGenerator gen(orbitals, userInput);

        // Generate integrals
        pHFOperator hf = basis_generator.GetClosedHFOperator();
        pHFIntegrals hf_electron(new HFIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

        pCoulombOperator coulomb(new CoulombOperator(lattice));
        pHartreeY hartreeY(new HartreeY(hf->GetIntegrator(), coulomb));
        pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron(new TwoElectronCoulombOperator<pSlaterIntegrals>(integrals));

        pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();

        pHamiltonianID key = std::make_shared<HamiltonianID>(sym);
        pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(sym, angular_library);
        HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
        H.GenerateMatrix();
        LevelVector levels = H.SolveMatrix(key, 3);

        // Calculate transition
        pIntegrator integrator(new SimpsonsIntegrator(lattice));
        pTimeDependentSpinorOperator M1 = std::make_shared<MJOperator>(1, integrator);
        M1->SetFrequency(levels[2]->GetEnergy() - levels[1]->GetEnergy());
        pTransitionIntegrals M1_matrix_elements(new TransitionIntegrals(orbitals, M1));
        M1_matrix_elements->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
        ManyBodyOperator<pTransitionIntegrals> M1_many_body(M1_matrix_elements);

        m1_hole = M1_many_body.GetMatrixElement(*levels[2], *levels[1]);

        // Convert to strength = (reduced matrix element)^2
        m1_hole = m1_hole/math->Electron3j(sym.GetTwoJ(), sym.GetTwoJ(), 1, sym.GetTwoJ(), -sym.GetTwoJ());
        m1_hole = m1_hole * m1_hole;
    }

    EXPECT_NEAR(m1_hole, m1_electron, 1.e-6 * fabs(m1_hole));
}

TEST(EJOperatorTester, Screening)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Ne
    std::string user_input_string = std::string() +
        "NuclearRadius = 3.5\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 10\n" +
        "[HF]\n" +
        "N = 10\n" +
        "Configuration = '1s2 2s2 2p6'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 3sp\n" +
        "BSpline/Rmax = 50.0\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();

    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pSpinorOperator pE1 = std::make_shared<EJOperator>(1, integrator);

    DebugOptions.LogHFIterations(true);
    pRPASolver rpa_solver = std::make_shared<RPASolver>(lattice);
    pRPAOperator rpa = std::make_shared<RPAOperator>(pE1, basis_generator.GetClosedHFOperator(), basis_generator.GetHartreeY(), rpa_solver);
    rpa->SetFrequency(0.0);

    pOrbital s = std::make_shared<Orbital>(-1);
    s->resize(1000);
    s->f = std::vector<double>(1000, 1.);
    s->dfdr = std::vector<double>(1000, 0.);
    s->g = s->dfdr;
    s->dgdr = s->dfdr;

    SpinorFunction fs = pE1->ApplyTo(*s, 1);
    SpinorFunction dVs = rpa->ApplyTo(*s, 1);

    // Expect Z_ion/Z: see, e.g., Dzuba et al. Phys. Lett. A 118, 177 (1986)
    EXPECT_NEAR(0., dVs.f[0]/fs.f[0], 1.e-4);
}
