#include "HamiltonianMatrix.h"
#include "gtest/gtest.h"
#include "Include.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Basis/BasisGenerator.h"
#include "ConfigGenerator.h"
#include "GFactor.h"
#include "Atom/MultirunOptions.h"

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
    pOneElectronIntegrals hf_electron(new OneElectronIntegrals(orbitals, hf));
    hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

    pCoulombOperator coulomb(new CoulombOperator(lattice));
    pHartreeY hartreeY(new HartreeY(hf->GetOPIntegrator(), coulomb));
    pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
    integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);

    EXPECT_EQ(31747, integrals->size());

    ConfigGenerator config_generator(orbitals, userInput);
    pRelativisticConfigList relconfigs;
    Symmetry sym(0, Parity::even);
    pLevelMap levels = pLevelMap(new LevelMap());
    double ground_state_energy = 0., excited_state_energy = 0.;

    // Generate matrix and configurations
    relconfigs = config_generator.GenerateRelativisticConfigurations(sym);
    pTwoElectronCoulombOperator twobody_electron(new TwoElectronCoulombOperator<pSlaterIntegrals>(integrals));
    HamiltonianMatrix H_even(hf_electron, twobody_electron, relconfigs);
    H_even.GenerateMatrix();

    // Solve matrix
    H_even.SolveMatrix(sym, 3, levels);
    levels->Print(sym);
    ground_state_energy = levels->begin(sym)->second->GetEnergy();

    // Rinse and repeat for excited state
    sym = Symmetry(2, Parity::odd);
    relconfigs = config_generator.GenerateRelativisticConfigurations(sym);

    HamiltonianMatrix H_odd(hf_electron, twobody_electron, relconfigs);
    H_odd.GenerateMatrix();
    H_odd.SolveMatrix(sym, 3, levels);
    GFactorCalculator g_factors(hf->GetOPIntegrator(), orbitals);
    g_factors.CalculateGFactors(*levels, sym);
    levels->Print(sym);
    excited_state_energy = levels->begin(sym)->second->GetEnergy();

    EXPECT_NEAR(1.5, levels->begin(sym)->second->GetgFactor(), 0.001);

    // Check energy 3s2 -> 3s3p J = 1 (should be within 20%)
    EXPECT_NEAR(21870, (excited_state_energy - ground_state_energy) * MathConstant::Instance()->HartreeEnergyInInvCm(), 4000);
}

TEST(HamiltonianMatrixTester, HolesOnly)
{
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);

    pLattice lattice(new Lattice(1000, 1.e-6, 50.));
    std::vector<Symmetry> symmetries;
    symmetries.emplace_back(2, Parity::even);
    symmetries.emplace_back(4, Parity::even);

    pLevelMap electron_levels(new LevelMap());
    pLevelMap hole_levels(new LevelMap());

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
        pOneElectronIntegrals hf_electron(new OneElectronIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

        pCoulombOperator coulomb(new CoulombOperator(lattice));
        pHartreeY hartreeY(new HartreeY(hf->GetOPIntegrator(), coulomb));
        pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron(new TwoElectronCoulombOperator<pSlaterIntegrals>(integrals));
        GFactorCalculator g_factors(hf->GetOPIntegrator(), orbitals);

        for(auto sym : symmetries)
        {
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(sym);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            H.SolveMatrix(sym, 6, electron_levels);
            g_factors.CalculateGFactors(*electron_levels, sym);
            electron_levels->Print(sym);
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
        pOneElectronIntegrals hf_electron(new OneElectronIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

        pCoulombOperator coulomb(new CoulombOperator(lattice));
        pHartreeY hartreeY(new HartreeY(hf->GetOPIntegrator(), coulomb));
        pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron(new TwoElectronCoulombOperator<pSlaterIntegrals>(integrals));
        GFactorCalculator g_factors(hf->GetOPIntegrator(), orbitals);

        for(auto sym : symmetries)
        {
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(sym);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            H.SolveMatrix(sym, 6, hole_levels);
            g_factors.CalculateGFactors(*hole_levels, sym);
            hole_levels->Print(sym);
        }
    }

    ASSERT_EQ(electron_levels->size(), hole_levels->size());

    LevelID base_level = LevelID(symmetries[0], 0);
    double electron_base = (*electron_levels)[base_level]->GetEnergy();
    double hole_base = (*hole_levels)[base_level]->GetEnergy();

    for(auto& pair: *electron_levels)
    {
        if(pair.first != base_level)
        {
            auto hole_it = hole_levels->find(pair.first);
            double e = pair.second->GetEnergy() - electron_base;
            double h = hole_it->second->GetEnergy() - hole_base;

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

    pLevelMap electron_levels(new LevelMap());
    pLevelMap hole_levels(new LevelMap());

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
        pOneElectronIntegrals hf_electron(new OneElectronIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

        pCoulombOperator coulomb(new CoulombOperator(lattice));
        pHartreeY hartreeY(new HartreeY(hf->GetOPIntegrator(), coulomb));
        pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron(new TwoElectronCoulombOperator<pSlaterIntegrals>(integrals));
        GFactorCalculator g_factors(hf->GetOPIntegrator(), orbitals);

        for(auto sym : symmetries)
        {
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(sym);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            H.SolveMatrix(sym, 6, electron_levels);
            g_factors.CalculateGFactors(*electron_levels, sym);
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
        pOneElectronIntegrals hf_electron(new OneElectronIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);

        pCoulombOperator coulomb(new CoulombOperator(lattice));
        pHartreeY hartreeY(new HartreeY(hf->GetOPIntegrator(), coulomb));
        pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron(new TwoElectronCoulombOperator<pSlaterIntegrals>(integrals));
        GFactorCalculator g_factors(hf->GetOPIntegrator(), orbitals);

        for(auto sym : symmetries)
        {
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(sym);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            H.SolveMatrix(sym, 6, hole_levels);
            g_factors.CalculateGFactors(*hole_levels, sym);
        }
    }

    ASSERT_EQ(electron_levels->size(), hole_levels->size());

    LevelID base_level = LevelID(symmetries[0], 0);
    double electron_base = (*electron_levels)[base_level]->GetEnergy();
    double hole_base = (*hole_levels)[base_level]->GetEnergy();

    for(auto& pair: *electron_levels)
    {
        if(pair.first != base_level)
        {
            auto hole_it = hole_levels->find(pair.first);
            double e = pair.second->GetEnergy() - electron_base;
            double h = hole_it->second->GetEnergy() - hole_base;

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
    
    pLevelMap electron_levels(new LevelMap());
    pLevelMap hole_levels(new LevelMap());
    
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
        pOneElectronIntegrals hf_electron(new OneElectronIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
        
        pCoulombOperator coulomb(new CoulombOperator(lattice));
        pHartreeY hartreeY(new HartreeY(hf->GetOPIntegrator(), coulomb));
        pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron(new TwoElectronCoulombOperator<pSlaterIntegrals>(integrals));
        GFactorCalculator g_factors(hf->GetOPIntegrator(), orbitals);
        
        for(auto sym : symmetries)
        {
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(sym);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            H.SolveMatrix(sym, 6, electron_levels);
            g_factors.CalculateGFactors(*electron_levels, sym);
            electron_levels->Print(sym);
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
        "FrozenCore = \n" +
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
        pOneElectronIntegrals hf_electron(new OneElectronIntegrals(orbitals, hf));
        hf_electron->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
        
        pCoulombOperator coulomb(new CoulombOperator(lattice));
        pHartreeY hartreeY(new HartreeY(hf->GetOPIntegrator(), coulomb));
        pSlaterIntegrals integrals(new SlaterIntegralsMap(orbitals, hartreeY));
        integrals->CalculateTwoElectronIntegrals(orbitals->valence, orbitals->valence, orbitals->valence, orbitals->valence);
        pTwoElectronCoulombOperator twobody_electron(new TwoElectronCoulombOperator<pSlaterIntegrals>(integrals));
        GFactorCalculator g_factors(hf->GetOPIntegrator(), orbitals);
        
        for(auto sym : symmetries)
        {
            pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(sym);
            HamiltonianMatrix H(hf_electron, twobody_electron, relconfigs);
            H.GenerateMatrix();
            H.SolveMatrix(sym, 6, hole_levels);
            g_factors.CalculateGFactors(*hole_levels, sym);
            hole_levels->Print(sym);
        }
    }
    
    ASSERT_EQ(electron_levels->size(), hole_levels->size());
    
    LevelID base_level = LevelID(symmetries[0], 0);
    double electron_base = (*electron_levels)[base_level]->GetEnergy();
    double hole_base = (*hole_levels)[base_level]->GetEnergy();
    
    for(auto& pair: *electron_levels)
    {
        if(pair.first != base_level)
        {
            auto hole_it = hole_levels->find(pair.first);
            double e = pair.second->GetEnergy() - electron_base;
            double h = hole_it->second->GetEnergy() - hole_base;
            
            EXPECT_NEAR(e, h, mmax(1.e-5 * fabs(e), 1.e-12));
        }
    }
}
