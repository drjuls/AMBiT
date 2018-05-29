#include "ManyBodyOperator.h"
#include "gtest/gtest.h"
#include "Include.h"

using namespace Ambit;

TEST(ManyBodyOperatorTester, ProjectionDifferencesElectrons)
{
    // Make Projections and check projection differences
    // Remember: Relativistic configurations are sorted first on holes, then abs(kappa), then kappa, then pqn

    ManyBodyOperator<> many_body_operator;
    RelativisticConfiguration config1, config2;
    std::vector<int> TwoM1, TwoM2;
    ManyBodyOperator<>::IndirectProjectionStruct indirects;
    std::array<std::pair<OrbitalInfo, OrbitalInfo>, 2> diffarray;
    std::array<std::pair<OrbitalInfo, OrbitalInfo>, 1> diffarray1;

    // < e1 e2 | G | e3 e4 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, -2), 1));
        config2.insert(std::make_pair(OrbitalInfo(3, -1), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, -3), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(-1); TwoM1.push_back(3);
        TwoM2.push_back(1); TwoM2.push_back(1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(2, diffs);

        ASSERT_EQ(abs(diffs), config1.GetConfigDifferences<2>(config2, diffarray));
        for(int i = 0; i < abs(diffs); i++)
        {   EXPECT_EQ(OrbitalInfo(*indirects.left[i]), diffarray[i].first);
            EXPECT_EQ(OrbitalInfo(*indirects.right[i]), diffarray[i].second);
        }
    }

    // < e1 e2 | G | e1 e3 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, -2), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, -3), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(-1); TwoM1.push_back(3);
        TwoM2.push_back(-1); TwoM2.push_back(3);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(1, diffs);

        ASSERT_EQ(abs(diffs), config1.GetConfigDifferences<1>(config2, diffarray1));
        for(int i = 0; i < abs(diffs); i++)
        {   EXPECT_EQ(OrbitalInfo(*indirects.left[i]), diffarray1[i].first);
            EXPECT_EQ(OrbitalInfo(*indirects.right[i]), diffarray1[i].second);
        }
    }

    // < e1 e2 | G | e3 e1 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, -2), 1));
        config2.insert(std::make_pair(OrbitalInfo(3, -1), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, 1), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(-1); TwoM1.push_back(3);
        TwoM2.push_back(3); TwoM2.push_back(-1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(-1, diffs);

        ASSERT_EQ(abs(diffs), config1.GetConfigDifferences<1>(config2, diffarray1));
        for(int i = 0; i < abs(diffs); i++)
        {   EXPECT_EQ(OrbitalInfo(*indirects.left[i]), diffarray1[i].first);
            EXPECT_EQ(OrbitalInfo(*indirects.right[i]), diffarray1[i].second);
        }

        many_body_operator.make_indirect_projection(proj1, indirects.right);
        many_body_operator.make_indirect_projection(proj2, indirects.left);
        diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(-1, diffs);
    }

    // < e1 e2 e3 | G | e4 e5 e1 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, -2), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, 2), 1));
        config2.insert(std::make_pair(OrbitalInfo(3, -1), 2));
        config2.insert(std::make_pair(OrbitalInfo(4, 1), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(-1); TwoM1.push_back(3); TwoM1.push_back(-3);
        TwoM2.push_back(-1); TwoM2.push_back(1); TwoM2.push_back(-1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        EXPECT_EQ(2, diffs);

        ASSERT_EQ(abs(diffs), config1.GetConfigDifferences<2>(config2, diffarray));
        for(int i = 0; i < abs(diffs); i++)
        {   EXPECT_EQ(OrbitalInfo(*indirects.left[i]), diffarray[i].first);
            EXPECT_EQ(OrbitalInfo(*indirects.right[i]), diffarray[i].second);
        }
    }

    // < e1 e2 e3 | G | e4 e1 e5 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, -2), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, 2), 1));
        config2.insert(std::make_pair(OrbitalInfo(3, -1), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, 2), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(-1); TwoM1.push_back(3); TwoM1.push_back(-3);
        TwoM2.push_back(-1); TwoM2.push_back(-1); TwoM2.push_back(1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        EXPECT_EQ(-2, diffs);

        ASSERT_EQ(abs(diffs), config1.GetConfigDifferences<2>(config2, diffarray));
        for(int i = 0; i < abs(diffs); i++)
        {   EXPECT_EQ(OrbitalInfo(*indirects.left[i]), diffarray[i].first);
            EXPECT_EQ(OrbitalInfo(*indirects.right[i]), diffarray[i].second);
        }
    }
}

TEST(ManyBodyOperatorTester, ProjectionDifferencesHoles)
{
    // Make Projections and check projection differences
    // Remember: Relativistic configurations are sorted first on hole, then abs(kappa), then kappa, then pqn

    ManyBodyOperator<> many_body_operator;
    RelativisticConfiguration config1, config2;
    std::vector<int> TwoM1, TwoM2;
    ManyBodyOperator<>::IndirectProjectionStruct indirects;
    std::array<std::pair<OrbitalInfo, OrbitalInfo>, 2> diffarray;
    std::array<std::pair<OrbitalInfo, OrbitalInfo>, 3> diffarray3;

    // < h3 e2 | G | h1 e4 > => - < h1 e2 | G | h3 e4 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(3, -2), -1));
        config1.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config2.insert(std::make_pair(OrbitalInfo(3, -1), -1));
        config2.insert(std::make_pair(OrbitalInfo(4, -3), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(-1); TwoM1.push_back(3);
        TwoM2.push_back(1); TwoM2.push_back(1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(-2, diffs);

        ASSERT_EQ(abs(diffs), config1.GetConfigDifferences<2>(config2, diffarray));
        for(int i = 0; i < abs(diffs); i++)
        {   EXPECT_EQ(OrbitalInfo(*indirects.left[i]), diffarray[i].first);
            EXPECT_EQ(OrbitalInfo(*indirects.right[i]), diffarray[i].second);
        }
    }

    // < e1 | G | h2 e2 e3> => - < h2 e1 | G | e2 e3 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, -1), 1));

        config2.insert(std::make_pair(OrbitalInfo(4, -3), -1));
        config2.insert(std::make_pair(OrbitalInfo(3, -2), 2));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(1);
        TwoM2.push_back(1); TwoM2.push_back(-1); TwoM2.push_back(1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(-2, diffs);

        ASSERT_EQ(abs(diffs), config1.GetConfigDifferences<2>(config2, diffarray));
        for(int i = 0; i < abs(diffs); i++)
        {   EXPECT_EQ(OrbitalInfo(*indirects.left[i]), diffarray[i].first);
            EXPECT_EQ(OrbitalInfo(*indirects.right[i]), diffarray[i].second);
        }
    }

    // < e1 | G | h2 e1 e3>
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, -1), 1));

        config2.insert(std::make_pair(OrbitalInfo(4, -1), 1));
        config2.insert(std::make_pair(OrbitalInfo(3, -2), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, -3), -1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(1);
        TwoM2.push_back(-1); TwoM2.push_back(1); TwoM2.push_back(1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(1, diffs);

        ASSERT_EQ(abs(diffs), config1.GetConfigDifferences<2>(config2, diffarray));
        for(int i = 0; i < abs(diffs); i++)
        {   EXPECT_EQ(OrbitalInfo(*indirects.left[i]), diffarray[i].first);
            EXPECT_EQ(OrbitalInfo(*indirects.right[i]), diffarray[i].second);
        }
    }

    // < h1 h2 e1 e2 | G | 0 > => - < e1 e2 | G | h1 h2 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, -1), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, -2), -1));
        config1.insert(std::make_pair(OrbitalInfo(4, -3), -1));
        config1.insert(std::make_pair(OrbitalInfo(4, -4), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(1); TwoM1.push_back(1); TwoM1.push_back(1); TwoM1.push_back(1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(-2, diffs);

        ASSERT_EQ(abs(diffs), config1.GetConfigDifferences<2>(config2, diffarray));
        for(int i = 0; i < abs(diffs); i++)
        {   EXPECT_EQ(OrbitalInfo(*indirects.left[i]), diffarray[i].first);
            EXPECT_EQ(OrbitalInfo(*indirects.right[i]), diffarray[i].second);
        }

        many_body_operator.make_indirect_projection(proj1, indirects.right);
        many_body_operator.make_indirect_projection(proj2, indirects.left);
        diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        EXPECT_EQ(-2, diffs);
    }

    // < h1 h2 e1 | G | h3 h4 e2 > => < h3 h4 e1 | G | h1 h2 e2>
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, -1), -1));
        config1.insert(std::make_pair(OrbitalInfo(3, -2), -1));
        config1.insert(std::make_pair(OrbitalInfo(5, -1), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, -3), -1));
        config2.insert(std::make_pair(OrbitalInfo(4, -4), -1));
        config2.insert(std::make_pair(OrbitalInfo(5, 1), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(1); TwoM1.push_back(1); TwoM1.push_back(1);
        TwoM2.push_back(1); TwoM2.push_back(1); TwoM2.push_back(1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<3>(indirects);
        EXPECT_EQ(3, diffs);

        ASSERT_EQ(abs(diffs), config1.GetConfigDifferences<3>(config2, diffarray3));
        for(int i = 0; i < abs(diffs); i++)
        {   EXPECT_EQ(OrbitalInfo(*indirects.left[i]), diffarray3[i].first);
            EXPECT_EQ(OrbitalInfo(*indirects.right[i]), diffarray3[i].second);
        }
    }
}

TEST(ManyBodyOperatorTester, ProjectionDifferencesSkipped)
{
    // Make Projections and check projection differences
    // Remember: Relativistic configurations are sorted first on hole, then abs(kappa), then kappa, then pqn
    // but skipped p2 (right projection) electron is last (arbitrarily, always).

    ManyBodyOperator<> many_body_operator;
    RelativisticConfiguration config1, config2;
    std::vector<int> TwoM1, TwoM2;
    ManyBodyOperator<>::IndirectProjectionStruct indirects;
    std::array<std::pair<OrbitalInfo, OrbitalInfo>, 2> diffarray;

    // Electrons only:
    // < e1 e2 | G | e3 ek > => - < e1 e2 | G | ek e3 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, -2), 1));
        config2.insert(std::make_pair(OrbitalInfo(3, -1), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(-1); TwoM1.push_back(3);
        TwoM2.push_back(1);

        ElectronInfo skipped(13, -3, 1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects, &skipped);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(-2, diffs);

        EXPECT_EQ(abs(diffs), config1.GetConfigDifferencesCount(config2));
    }

    // < e1 e2 | G | e1 ek >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, -2), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, 1), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(-1); TwoM1.push_back(3);
        TwoM2.push_back(-1);

        ElectronInfo skipped(13, -3, 3);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects, &skipped);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(1, diffs);

        EXPECT_EQ(abs(diffs), config1.GetConfigDifferencesCount(config2));
    }

    // < e1 e2 e3 | G | e4 e1 ek > => < e2 e3 | g | ek e4 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, -2), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, 2), 1));
        config2.insert(std::make_pair(OrbitalInfo(3, -1), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, 1), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(-1); TwoM1.push_back(3); TwoM1.push_back(-3);
        TwoM2.push_back(-1); TwoM2.push_back(-1);

        ElectronInfo skipped(13, -1, 1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects, &skipped);
        EXPECT_EQ(2, diffs);

        EXPECT_EQ(abs(diffs), config1.GetConfigDifferencesCount(config2));
    }

    // < h3 e2 | G | h1 ek > => - < h1 e2 | G | h3 ek > = < h1 e2 | G | ek h3 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(3, -2), -1));
        config1.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config2.insert(std::make_pair(OrbitalInfo(3, -1), -1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(-1); TwoM1.push_back(3);
        TwoM2.push_back(1);

        ElectronInfo skipped(13, -3, 1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects, &skipped);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(2, diffs);

        EXPECT_EQ(abs(diffs), config1.GetConfigDifferencesCount(config2));
    }

    // < e1 | G | h2 e2 ek> => - < h2 e1 | G | e2 ek > = < h2 e1 | G | ek e2 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, -1), 1));

        config2.insert(std::make_pair(OrbitalInfo(4, -3), -1));
        config2.insert(std::make_pair(OrbitalInfo(3, -2), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(1);
        TwoM2.push_back(1); TwoM2.push_back(-1);

        ElectronInfo skipped(13, -2, 1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects, &skipped);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(2, diffs);

        EXPECT_EQ(abs(diffs), config1.GetConfigDifferencesCount(config2));
    }

    // < e1 | G | h2 e1 ek> => < h2 e1 | G | ek e1 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, -1), 1));

        config2.insert(std::make_pair(OrbitalInfo(4, -3), -1));
        config2.insert(std::make_pair(OrbitalInfo(4, -1), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(1);
        TwoM2.push_back(-1); TwoM2.push_back(1);

        ElectronInfo skipped(13, -2, 1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects, &skipped);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(1, diffs);

        EXPECT_EQ(abs(diffs), config1.GetConfigDifferencesCount(config2));
    }

    // < 0 | G | h1 h2 e1 ek > => - < h1 h2 | G | e1 ek > = < h1 h2 | G | ek e1 >
    {   config1.clear(); config2.clear();
        config2.insert(std::make_pair(OrbitalInfo(4, -1), 1));
        config2.insert(std::make_pair(OrbitalInfo(3, -2), -1));
        config2.insert(std::make_pair(OrbitalInfo(4, -3), -1));

        TwoM1.clear(); TwoM2.clear();
        TwoM2.push_back(1); TwoM2.push_back(1); TwoM2.push_back(1);

        ElectronInfo skipped(11, -4, 1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects, &skipped);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(2, diffs);

        EXPECT_EQ(abs(diffs), config1.GetConfigDifferencesCount(config2));
    }

    // < h1 h2 e1 | G | h3 h1 ek > => < h3 e1 | G | h2 ek > =  - < h3 e1 | G | ek h2 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, -2), -1));
        config1.insert(std::make_pair(OrbitalInfo(3, -3), -1));
        config1.insert(std::make_pair(OrbitalInfo(5, -1), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, -1), -1));
        config2.insert(std::make_pair(OrbitalInfo(4, -2), -1));

        TwoM1.clear(); TwoM2.clear();
        TwoM1.push_back(1); TwoM1.push_back(1); TwoM1.push_back(1);
        TwoM2.push_back(1); TwoM2.push_back(1);

        ElectronInfo skipped(13, -1, 1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects, &skipped);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(-2, diffs);

        EXPECT_EQ(abs(diffs), config1.GetConfigDifferencesCount(config2));
    }
}
