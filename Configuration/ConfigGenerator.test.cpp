#include "ConfigGenerator.h"
#include "gtest/gtest.h"
#include "Include.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/ConfigurationParser.h"

TEST(ConfigGeneratorTester, CountConfigurations)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // CaII
    OccupationMap occ = ConfigurationParser::ParseFractionalConfiguration("1s2 2s2 2p6 3s2 3p6");

    MultirunOptions userInput("template.input", "//", "\n", ",");

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

    int total_non_rel = 0;
    int total_rel = 0;
    ConfigGenerator gen(core, excited, userInput);

    pRelativisticConfigList relconfigs = gen.GenerateRelativisticConfigurations(Symmetry(0, Parity::even));
    ConfigList nonrelconfigs(*relconfigs);
    total_non_rel += nonrelconfigs.size();
    total_rel += relconfigs->size();

    relconfigs = gen.GenerateRelativisticConfigurations(Symmetry(0, Parity::odd));
    nonrelconfigs = *relconfigs;
    total_non_rel += nonrelconfigs.size();
    total_rel += relconfigs->size();

    // Check number of non-rel configurations:
    // 2 electrons into 4-6s, 4-6p, 3-6d, 4-6f: total 13 non-rel states
    //  = (13 choose 2) + 13 (doubles) = 91
    EXPECT_EQ(91, total_non_rel);

    // Number of relativistic configurations:
    // 2 electrons into 23 orbitals
    EXPECT_EQ(276, total_rel);

    // Test RelativisticConfigList projection iterator
    int projection_count = 0;
    for(auto rconfig: *relconfigs)
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
