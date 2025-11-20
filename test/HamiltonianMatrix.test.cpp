#include "Configuration/HamiltonianMatrix.h"
#include "gtest/gtest.h"
#include <memory>
#include "Include.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Basis/BasisGenerator.h"
#include "Configuration/ConfigGenerator.h"
#include "Configuration/GFactor.h"
#include "Configuration/LevelMap.h"
#include "Atom/MultirunOptions.h"

using namespace Ambit;

TEST(HamiltonianMatrixTester, MgILevels)
{
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);

    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // MgI
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
        "[CI]\n" +
        "LeadingConfigurations = '3s2, 3s1 3p1'\n" +
        "ElectronExcitations = 1\n" +
        "EvenParityTwoJ = '0'\n" +
        "OddParityTwoJ = '0, 2'\n" +
        "NumSolutions = 3\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();

    // Generate integrals
    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pHFIntegrals hf_electron(new HFIntegrals(orbitals, hf));
    hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

    pCoulombOperator coulomb(new CoulombOperator(lattice));
    pHartreeY hartreeY(new HartreeY(hf->GetIntegrator(), coulomb));
    pSlaterIntegrals integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
    integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);

    EXPECT_EQ(31747, integrals->size());

    ConfigGenerator config_generator(orbitals, userInput);
    pRelativisticConfigList relconfigs;
    Symmetry sym(0, Parity::even);
    double ground_state_energy = 0., excited_state_energy = 0.;

    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();

    // Generate matrix and configurations
    auto configs = config_generator.GenerateConfigurations();
    relconfigs = config_generator.GenerateRelativisticConfigurations(configs, sym, angular_library);
    pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
    HamiltonianMatrix H_even(hf_electron, twobody_electron, relconfigs);
    H_even.GenerateMatrix();

    // Solve matrix
    LevelVector levels = H_even.SolveMatrix(std::make_shared<HamiltonianID>(sym), 3);
    levels.Print();
    ground_state_energy = levels.eigenvalues[0];

    // Rinse and repeat for excited state
    sym = Symmetry(2, Parity::odd);
    relconfigs = config_generator.GenerateRelativisticConfigurations(configs, sym, angular_library);

    HamiltonianMatrix H_odd(hf_electron, twobody_electron, relconfigs);
    H_odd.GenerateMatrix();
    levels = H_odd.SolveMatrix(std::make_shared<HamiltonianID>(sym), 3);
    GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);
    g_factors.CalculateGFactors(levels);
    levels.Print();
    excited_state_energy = levels.eigenvalues[0];

    EXPECT_NEAR(1.5, levels.g_factors[0], 0.001);

    // Check energy 3s2 -> 3s3p J = 1 (should be within 20%)
    EXPECT_NEAR(21870, (excited_state_energy - ground_state_energy) * MathConstant::Instance()->HartreeEnergyInInvCm(), 4000);
}

TEST(HamiltonianMatrixTester, MgISaveMatrix)
{
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);

    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // MgI
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
        "[CI]\n" +
        "LeadingConfigurations = '3s2, 3s1 3p1'\n" +
        "ElectronExcitations = 1\n" +
        "EvenParityTwoJ = '0'\n" +
        "OddParityTwoJ = '0, 2'\n" +
        "NumSolutions = 3\n" +
        "Output/--write-hamiltonian\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    std::string even_filename = "MgITest.e.matrix";
    std::string odd_filename = "MgITest.o.matrix";

    // Energies we got when (w)riting and (r)eading. These should be identical, as reading the
    // matrix from the file should have the same results as generating from scratch
    double wground_state_energy = 0., wexcited_state_energy = 0.; // Energy when writing
    double rground_state_energy = 0., rexcited_state_energy = 0.; // Energy when reading

    // Generate matrix and configurations, then save to disk
    {
        // Get core and excited basis
        BasisGenerator basis_generator(lattice, userInput);
        pCore core = basis_generator.GenerateHFCore();
        pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();

        // Generate integrals
        pHFOperator hf = basis_generator.GetClosedHFOperator();
        pHFIntegrals hf_electron(new HFIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

        pCoulombOperator coulomb(new CoulombOperator(lattice));
        pHartreeY hartreeY(new HartreeY(hf->GetIntegrator(), coulomb));
        pSlaterIntegrals integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        ConfigGenerator config_generator(orbitals, userInput);

        pRelativisticConfigList relconfigs;
        Symmetry sym(0, Parity::even);

        pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();

        auto configs = config_generator.GenerateConfigurations();
        relconfigs = config_generator.GenerateRelativisticConfigurations(configs, sym, angular_library);
        pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
        auto H_even = std::make_shared<HamiltonianMatrix>(hf_electron, twobody_electron, relconfigs);
        H_even->GenerateMatrix();
        // Write to disk
        H_even->Write(even_filename);
        // Solve matrix
        LevelVector levels = H_even->SolveMatrix(std::make_shared<HamiltonianID>(sym), 3);
        levels.Print();
        wground_state_energy = levels.eigenvalues[0];

        // Rinse and repeat for excited state
        sym = Symmetry(2, Parity::odd);
        relconfigs = config_generator.GenerateRelativisticConfigurations(configs, sym, angular_library);

        auto H_odd = std::make_shared<HamiltonianMatrix> (hf_electron, twobody_electron, relconfigs);
        H_odd->GenerateMatrix();
        // Write to disk
        H_odd->Write(odd_filename);
        // Solve
        levels = H_odd->SolveMatrix(std::make_shared<HamiltonianID>(sym), 3);
        GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);
        g_factors.CalculateGFactors(levels);
        levels.Print();
        wexcited_state_energy = levels.eigenvalues[0];
    }
    // Now read back from disk and check we get the correct values
    {
        // Get core and excited basis
        BasisGenerator basis_generator(lattice, userInput);
        pCore core = basis_generator.GenerateHFCore();
        pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();

        // Generate integrals
        pHFOperator hf = basis_generator.GetClosedHFOperator();
        pHFIntegrals hf_electron(new HFIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

        pCoulombOperator coulomb(new CoulombOperator(lattice));
        pHartreeY hartreeY(new HartreeY(hf->GetIntegrator(), coulomb));
        pSlaterIntegrals integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        ConfigGenerator config_generator(orbitals, userInput);

        pRelativisticConfigList relconfigs;
        Symmetry sym(0, Parity::even);

        pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();

        auto configs = config_generator.GenerateConfigurations();
        relconfigs = config_generator.GenerateRelativisticConfigurations(configs, sym, angular_library);
        pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
        auto H_even = std::make_shared<HamiltonianMatrix>(hf_electron, twobody_electron, relconfigs);
        bool read_success;
        read_success = H_even->Read(even_filename);
        EXPECT_TRUE(read_success);
        // Solve matrix
        LevelVector levels = H_even->SolveMatrix(std::make_shared<HamiltonianID>(sym), 3);
        levels.Print();
        rground_state_energy = levels.eigenvalues[0];

        // Rinse and repeat for excited state
        sym = Symmetry(2, Parity::odd);
        relconfigs = config_generator.GenerateRelativisticConfigurations(configs, sym, angular_library);

        auto H_odd = std::make_shared<HamiltonianMatrix>(hf_electron, twobody_electron, relconfigs);
        read_success = H_odd->Read(odd_filename);
        EXPECT_TRUE(read_success);
        // Solve
        levels = H_odd->SolveMatrix(std::make_shared<HamiltonianID>(sym), 3);
        GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);
        g_factors.CalculateGFactors(levels);
        levels.Print();
        rexcited_state_energy = levels.eigenvalues[0];

        EXPECT_NEAR(1.5, levels.g_factors[0], 0.001);

        // Check energy 3s2 -> 3s3p J = 1 (should be within 20%)
        EXPECT_NEAR((wexcited_state_energy - wground_state_energy) * MathConstant::Instance()->HartreeEnergyInInvCm(), (rexcited_state_energy - rground_state_energy) * MathConstant::Instance()->HartreeEnergyInInvCm(), 4000);
    }
}
TEST(HamiltonianMatrixTester, HolesOnly)
{
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);

    pLattice lattice(new Lattice(1000, 1.e-6, 50.));
    std::vector<Symmetry> symmetries;
    symmetries.emplace_back(2, Parity::even);
    symmetries.emplace_back(4, Parity::even);

    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
    pLevelStore electron_levels = std::make_shared<LevelMap>(angular_library);
    pLevelStore hole_levels = std::make_shared<LevelMap>(angular_library);

    // CuIV - old electron way
    {   // Block for reusing names
        std::string user_input_string = std::string() +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 29\n" +
        "[HF]\n" +
        "N = 28\n" +
        "Configuration = '1s2 2s2 2p6: 3s2 3p6 3d10'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 3spd\n" +
        "BSpline/Rmax = 50.0\n" +
        "[CI]\n" +
        "LeadingConfigurations = '3s2 3p6 3d8'\n" +
        "ElectronExcitations = 2\n";

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
        pSlaterIntegrals integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
        GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);

        auto configs = gen.GenerateConfigurations();

        for(auto sym : symmetries)
        {
            pHamiltonianID key = std::make_shared<HamiltonianID>(sym);
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(configs, sym, angular_library);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            LevelVector levels = H.SolveMatrix(key, 6);
            g_factors.CalculateGFactors(levels);
            electron_levels->Store(key, levels);
            levels.Print();
        }
    }

    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(false);

    // CuIV - using holes now
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
        "ValenceBasis = 3spd\n" +
        "FrozenCore = 2sp\n" +
        "BSpline/Rmax = 50.0\n" +
        "[CI]\n" +
        "LeadingConfigurations = '3d-2'\n" +
        "ElectronExcitations = 0\n" +
        "HoleExcitations = 2\n";

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
        pSlaterIntegrals integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
        GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);

        auto configs = gen.GenerateConfigurations();

        for(auto sym : symmetries)
        {
            pHamiltonianID key = std::make_shared<HamiltonianID>(sym);
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(configs, sym, angular_library);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            LevelVector levels = H.SolveMatrix(key, 6);
            g_factors.CalculateGFactors(levels);
            hole_levels->Store(key, levels);
            levels.Print();
        }
    }

    ASSERT_EQ(electron_levels->keys.size(), hole_levels->keys.size());

    pHamiltonianID key = std::make_shared<HamiltonianID>(symmetries[0]);
    double electron_base = electron_levels->GetLevels(key).eigenvalues[0];
    double hole_base = hole_levels->GetLevels(key).eigenvalues[0];

    for(auto& key: *electron_levels)
    {
        LevelVector elv = electron_levels->GetLevels(key);
        LevelVector hlv = hole_levels->GetLevels(key);

        ASSERT_EQ(elv.NumLevels(), hlv.NumLevels());

        for(int i = 0; i < elv.NumLevels(); i++)
        {
            double e = elv.eigenvalues[i] - electron_base;
            double h = hlv.eigenvalues[i] - hole_base;

            EXPECT_NEAR(e, h, mmax(1.e-5 * fabs(e), 1.e-12));
        }
    }
}


TEST(HamiltonianMatrixTester, HolesVsElectrons)
{
    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(false);

    pLattice lattice(new Lattice(1000, 1.e-6, 50.));
    std::vector<Symmetry> symmetries;
    symmetries.emplace_back(1, Parity::even);
    symmetries.emplace_back(3, Parity::even);

    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
    pLevelStore electron_levels = std::make_shared<LevelMap>(angular_library);
    pLevelStore hole_levels = std::make_shared<LevelMap>(angular_library);

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
            "ValenceBasis = 4s3pd\n" +
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
        pSlaterIntegrals integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
        GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);

        pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
        auto configs = gen.GenerateConfigurations();

        for(auto sym : symmetries)
        {
            pHamiltonianID key = std::make_shared<HamiltonianID>(sym);
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(configs, sym, angular_library);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            LevelVector levels = H.SolveMatrix(key, 6);
            g_factors.CalculateGFactors(levels);
            electron_levels->Store(key, levels);
        }
    }

    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(false);

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
            "ValenceBasis = 4s3pd\n" +
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
        pSlaterIntegrals integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
        GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);

        pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
        auto configs = gen.GenerateConfigurations();

        for(auto sym : symmetries)
        {
            pHamiltonianID key = std::make_shared<HamiltonianID>(sym);
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(configs, sym, angular_library);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            LevelVector levels = H.SolveMatrix(key, 6);
            g_factors.CalculateGFactors(levels);
            hole_levels->Store(key, levels);
        }
    }

    ASSERT_EQ(electron_levels->keys.size(), hole_levels->keys.size());

    pHamiltonianID key = std::make_shared<HamiltonianID>(symmetries[0]);
    double electron_base = electron_levels->GetLevels(key).eigenvalues[0];
    double hole_base = hole_levels->GetLevels(key).eigenvalues[0];

    for(auto& key: *electron_levels)
    {
        LevelVector elv = electron_levels->GetLevels(key);
        LevelVector hlv = hole_levels->GetLevels(key);

        ASSERT_EQ(elv.NumLevels(), hlv.NumLevels());

        for(int i = 0; i < elv.NumLevels(); i++)
        {
            double e = elv.eigenvalues[i] - electron_base;
            double h = hlv.eigenvalues[i] - hole_base;

            EXPECT_NEAR(e, h, mmax(1.e-5 * fabs(e), 1.e-12));
        }
    }
}

TEST(HamiltonianMatrixTester, LiPlus)
{
    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(true);

    pLattice lattice(new Lattice(1000, 1.e-6, 50.));
    std::vector<Symmetry> symmetries;
    symmetries.emplace_back(0, Parity::even);

    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
    pLevelStore electron_levels = std::make_shared<LevelMap>(angular_library);
    pLevelStore hole_levels = std::make_shared<LevelMap>(angular_library);

    // LiII - old electron way
    {   // Block for reusing names
        std::string user_input_string = std::string() +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 3\n" +
        "[HF]\n" +
        "N = 2\n" +
        "Configuration = ':1s2'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 2s\n" +
        "BSpline/Rmax = 50.0\n" +
        "[CI]\n" +
        "LeadingConfigurations = '1s2'\n" +
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
        pSlaterIntegrals integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
        GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);

        pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
        auto configs = gen.GenerateConfigurations();

        for(auto sym : symmetries)
        {
            pHamiltonianID key = std::make_shared<HamiltonianID>(sym);
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(configs, sym, angular_library);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            LevelVector levels = H.SolveMatrix(key, 6);
            g_factors.CalculateGFactors(levels);
            electron_levels->Store(key, levels);
            levels.Print();
        }
    }

    // LiII - using holes now
    {   // Block for reusing names
        std::string user_input_string = std::string() +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 3\n" +
        "[HF]\n" +
        "N = 2\n" +
        "Configuration = '1s2'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 2s\n" +
        "FrozenCore = 0\n" +
        "BSpline/Rmax = 50.0\n" +
        "[CI]\n" +
        "LeadingConfigurations = '0'\n" +
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
        pSlaterIntegrals integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);
        GFactorCalculator g_factors(hf->GetIntegrator(), orbitals);

        pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
        auto configs = gen.GenerateConfigurations();

        for(auto sym : symmetries)
        {
            pHamiltonianID key = std::make_shared<HamiltonianID>(sym);
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(configs, sym, angular_library);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            LevelVector levels = H.SolveMatrix(key, 6);
            g_factors.CalculateGFactors(levels);
            hole_levels->Store(key, levels);
            levels.Print();
        }
    }

    ASSERT_EQ(electron_levels->keys.size(), hole_levels->keys.size());

    pHamiltonianID key = std::make_shared<HamiltonianID>(symmetries[0]);
    double electron_base = electron_levels->GetLevels(key).eigenvalues[0];
    double hole_base = hole_levels->GetLevels(key).eigenvalues[0];

    for(auto& key: *electron_levels)
    {
        LevelVector elv = electron_levels->GetLevels(key);
        LevelVector hlv = hole_levels->GetLevels(key);

        ASSERT_EQ(elv.NumLevels(), hlv.NumLevels());

        for(int i = 0; i < elv.NumLevels(); i++)
        {
            double e = elv.eigenvalues[i] - electron_base;
            double h = hlv.eigenvalues[i] - hole_base;

            EXPECT_NEAR(e, h, mmax(1.e-5 * fabs(e), 1.e-12));
        }
    }
}

TEST(HamiltonianMatrixTester, NonStretchedStates)
{
    // Testing levels when M != J
    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(false);

    pLattice lattice(new Lattice(1000, 1.e-6, 50.));
    std::vector<Symmetry> symmetries;
    symmetries.emplace_back(1, Parity::even);
    symmetries.emplace_back(3, Parity::even);

    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
    pLevelStore stretched_levels = std::make_shared<LevelMap>(angular_library);
    pLevelStore m_levels = std::make_shared<LevelMap>(angular_library);

    // CuIII - stretched state
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
        "ValenceBasis = 4spd\n" +
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
        pSlaterIntegrals integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);

        auto configs = gen.GenerateConfigurations();

        for(auto sym : symmetries)
        {
            pHamiltonianID key = std::make_shared<HamiltonianID>(sym);
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(configs, sym, angular_library);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            LevelVector levels = H.SolveMatrix(key, 6);
            stretched_levels->Store(key, levels);
        }
    }

    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(false);

    // CuIII - non-stretched state
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
        "ValenceBasis = 4spd\n" +
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
        pSlaterIntegrals integrals(new SlaterIntegralsFlatHash(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron = std::make_shared<TwoElectronCoulombOperator>(integrals);

        auto configs = gen.GenerateConfigurations();

        for(auto sym : symmetries)
        {
            pHamiltonianID key = std::make_shared<HamiltonianID>(sym);
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(configs, sym, angular_library);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            LevelVector levels = H.SolveMatrix(key, 6);
            m_levels->Store(key, levels);
        }
    }

    ASSERT_EQ(stretched_levels->keys.size(), m_levels->keys.size());

    for(auto& key: *stretched_levels)
    {
        LevelVector slv = stretched_levels->GetLevels(key);
        LevelVector mlv = m_levels->GetLevels(key);

        for(int i = 0; i < slv.NumLevels(); i++)
        {
            double e = slv.eigenvalues[i];
            double h = mlv.eigenvalues[i];

            EXPECT_NEAR(e, h, mmax(1.e-5 * fabs(e), 1.e-12));
        }
    }
}
