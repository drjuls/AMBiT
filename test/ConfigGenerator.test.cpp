#include "Configuration/ConfigGenerator.h"
#include "gtest/gtest.h"
#include "Include.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Basis/BasisGenerator.h"
#include "Specification/Specification.h"

using namespace Ambit;

TEST(ConfigGeneratorTester, CountConfigurations)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // CaII
    OccupationMap occ = ConfigurationParser::ParseFractionalConfiguration("1s2 2s2 2p6 3s2 3p6");

    std::string user_input_string = std::string() +
      "-c \n" +
      "ID = test\n" +
      "NuclearRadius = 3.7188\n" +
      "NuclearThickness = 2.3\n" +
      "Z = 20\n" +
      "[Lattice]\n" +
      "NumPoints = 1000\n" +
      "StartPoint = 1.e-6\n" +
      "EndPoint = 50.\n" +
      "[HF]\n" +
      "N = 18\n" +
      "Configuration = '1s2 2s2 2p6 3s2 3p6'\n" +
      "[Basis]\n" +
      "ValenceBasis = 6spdf\n" +
      "FrozenCore = 3s2p\n" +
      "BSpline/Rmax = 50.0\n" +
      "[CI]\n" +
      "LeadingConfigurations = '4s2, 4s1 4p1'\n" +
      "ExtraConfigurations   = '4p2'\n" +
      "ElectronExcitations = 2\n" +
      "HoleExcitations = 1\n" +
      "EvenParityTwoJ = '0'\n" +
      "OddParityTwoJ =  '0, 2'\n" +
      "NumSolutions = 3";

    auto specification = GlobalSpecification::Instance();
    // Can parse the user_input_string directly via parapara
    std::string perr = importSpecificationINIString(specification, user_input_string);
    if (!perr.empty()) {
        *errstream << "importSpecificationINIString:\n" << perr << std::endl;
        exit(1);
    }
    perr = specification->validateAndNormaliseSpecification();
    if (!perr.empty()) {
        *errstream << "validateAndNormaliseSpecification:\n" << perr << std::endl;
        exit(1);
    }

    pCore core(new Core(lattice));
    core->SetOccupancies(occ);

    pOrbitalMap excited = pOrbitalMap(new OrbitalMap(lattice));
    std::vector<int> limits = ConfigurationParser::ParseBasisSize("6spdf");
    for(int L = 0; L < limits.size(); L++)
    {
        for(int pqn = L + 1; pqn <= limits[L]; pqn++)
        {
            NonRelInfo info(pqn, L);

            OrbitalInfo orb_info_1 = info.GetFirstRelativisticInfo();
            if(core->GetState(orb_info_1) == NULL)
                excited->AddState(pOrbital(new Orbital(orb_info_1)));

            if(L > 0)
            {   OrbitalInfo orb_info_2 = info.GetSecondRelativisticInfo();
                if(core->GetState(orb_info_2) == NULL)
                    excited->AddState(pOrbital(new Orbital(orb_info_2)));
            }
        }
    }

    pOrbitalManager orbitals(new OrbitalManager(core, excited));

    int total_non_rel = 0;
    int total_rel = 0;
    ConfigGenerator gen(orbitals);

    pRelativisticConfigList configs = gen.GenerateConfigurations();
    total_non_rel += gen.GenerateNonRelConfigurations(configs)->first.size();

    pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(configs, Symmetry(0, Parity::even));
    total_rel += relconfigs->size();

    relconfigs = gen.GenerateRelativisticConfigurations(configs, Symmetry(0, Parity::odd));
    total_rel += relconfigs->size();

    // Check number of non-rel configurations:
    // 2 electrons into 4-6s, 4-6p, 3-6d, 4-6f: total 13 non-rel states
    //  = (13 choose 2) + 13 (doubles) = 91
    EXPECT_EQ(91, total_non_rel);

    // Number of relativistic configurations:
    // 2 electrons into 23 orbitals
    EXPECT_EQ(276, total_rel);

    // Test RelativisticConfigList projection iterator
    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
    gen.GenerateProjections(relconfigs, 0, angular_library);
    int projection_count = 0;
    for(auto& rconfig: *relconfigs)
    {
        projection_count += rconfig.projection_size();
    }

    int iterator_count = 0;
    auto proj_it = relconfigs->projection_begin();
    while(proj_it != relconfigs->projection_end())
    {
        iterator_count++;
        proj_it++;
    }

    EXPECT_EQ(projection_count, iterator_count);
    EXPECT_EQ(projection_count, relconfigs->projection_size());
}

TEST(ConfigGeneratorTester, HolesVsElectrons)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // CuII - old electron way
    std::string user_input_string = std::string() +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 29\n" +
        "[HF]\n" +
        "N = 28\n" +
        "Configuration = '1s2 2s2 2p6: 3s2 3p6 3d10'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 4spd\n" +
        "BSpline/Rmax = 50.0\n" +
        "[CI]\n" +
        "LeadingConfigurations = '3s2 3p6 3d8'\n" +
        "ElectronExcitations = 2\n" +
        "NumSolutions = 1\n";

    auto specification = GlobalSpecification::Instance();
    // Can parse the user_input_string directly via parapara
    std::string perr = importSpecificationINIString(specification, user_input_string);
    if (!perr.empty()) {
        *errstream << "importSpecificationINIString:\n" << perr << std::endl;
        exit(1);
    }

    perr = specification->validateAndNormaliseSpecification();
    if (!perr.empty()) {
        *errstream << "validateAndNormaliseSpecification:\n" << perr << std::endl;
        exit(1);
    }

    // Get core and excited basis
    BasisGenerator basis_generator(lattice);
    basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    ConfigGenerator gen(orbitals);

    auto configs = gen.GenerateConfigurations();
    pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(configs, Symmetry(0, Parity::even));
    unsigned int non_rel = gen.GenerateNonRelConfigurations(configs)->first.size();
    unsigned int rel = relconfigs->size();

    // CuII - using holes now
    std::string holes_input_string = std::string() +
        "NuclearRadius = 3.7188\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 29\n" +
        "[HF]\n" +
        "N = 28\n" +
        "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10'\n" +
        "[Basis]\n" +
        "--bspline-basis\n" +
        "ValenceBasis = 4spd\n" +
        "FrozenCore = 2sp\n" +
        "BSpline/Rmax = 50.0\n" +
        "[CI]\n" +
        "LeadingConfigurations = '3d-2'\n" +
        "ElectronExcitations = 2\n" +
        "HoleExcitations = 2\n" +
        "NumSolutions = 1\n";

    // Note: this will override the previsous specification, but we've already saved the values for
    // comparison
    // auto specification = GlobalSpecification::Instance();
    // Can parse the user_input_string directly via parapara
    specification->Reset();
    perr = importSpecificationINIString(specification, holes_input_string);
    if (!perr.empty()) {
        *errstream << "importSpecificationINIString:\n" << perr << std::endl;
        exit(1);
    }
    perr = specification->validateAndNormaliseSpecification();
    if (!perr.empty()) {
        *errstream << "validateAndNormaliseSpecification:\n" << perr << std::endl;
        exit(1);
    }

    // Get core and excited basis
    BasisGenerator holes_basis_generator(lattice);
    holes_basis_generator.GenerateHFCore();
    pOrbitalManagerConst holes_orbitals = holes_basis_generator.GenerateBasis();
    ConfigGenerator holes_gen(holes_orbitals);

    auto holes_configs = holes_gen.GenerateConfigurations();
    pRelativisticConfigList holes_relconfigs = holes_gen.GenerateRelativisticConfigurations(holes_configs, Symmetry(0, Parity::even));

    EXPECT_EQ(non_rel, holes_gen.GenerateNonRelConfigurations(holes_configs)->first.size());
    EXPECT_EQ(rel, holes_relconfigs->size());
}

TEST(ConfigGeneratorTester, NonSquare)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // CuIII
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
        "FrozenCore = 2sp\n" +
        "BSpline/Rmax = 50.0\n" +
        "[CI]\n" +
        "LeadingConfigurations = '3d-2'\n" +
        "ElectronExcitations = 2\n" +
        "HoleExcitations = 2\n" +
        "NumSolutions = 1\n" +
        "[CI/SmallSide]\n" +
        "LeadingConfigurations = '3d-2'\n" +
        "ElectronExcitations = 1 \n" +
        "ElectronBounds = '4spd' \n" +
        "HoleExcitations = 1\n" +
        "HoleBounds = '2sp'\n";

    auto specification = GlobalSpecification::Instance();
    // Can parse the user_input_string directly via parapara
    std::string perr = importSpecificationINIString(specification, user_input_string);
    if (!perr.empty()) {
        *errstream << "importSpecificationINIString:\n" << perr << std::endl;
        exit(1);
    }

    // Get core and excited basis
    BasisGenerator basis_generator(lattice);
    basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    ConfigGenerator gen(orbitals);

    auto configs = gen.GenerateConfigurations();
    pAngularDataLibrary angular_library = std::make_shared<AngularDataLibrary>();
    pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(configs, Symmetry(0, Parity::even), angular_library);
    unsigned int N = relconfigs->NumCSFs();
    unsigned int Nsmall = relconfigs->NumCSFsSmall();

    EXPECT_EQ(5322, N);
    EXPECT_EQ(43, Nsmall);
}
