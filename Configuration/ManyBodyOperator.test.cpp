#include "ManyBodyOperator.h"
#include "gtest/gtest.h"
#include "Include.h"

TEST(ManyBodyOperatorTester, ProjectionDifferencesPhase)
{
    // Make Projections and check projection differences
    // Remember: Relativistic configurations are sorted first on abs(kappa), then kappa, then pqn

    ManyBodyOperator<> many_body_operator;
    RelativisticConfiguration config1, config2;
    std::vector<int> TwoM1, TwoM2;
    ManyBodyOperator<>::IndirectProjectionStruct indirects;

    // < h3 e2 | G | h1 e4 >
    {   config1.clear(); config2.clear();
        config1.insert(std::make_pair(OrbitalInfo(4, 1), 1));
        config1.insert(std::make_pair(OrbitalInfo(3, -2), -1));
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
        EXPECT_EQ(2, diffs);
    }

    // < e1 | G | e1 h2 e3>
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
        EXPECT_EQ(2, diffs);
    }

    // < e1 h1 h2 e2 | G | 0 >
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
    }

    // < 0 | G | h1 e1 h2 e2 >
    {   config1.clear(); config2.clear();
        config2.insert(std::make_pair(OrbitalInfo(4, -1), -1));
        config2.insert(std::make_pair(OrbitalInfo(3, -2), 1));
        config2.insert(std::make_pair(OrbitalInfo(4, -3), -1));
        config2.insert(std::make_pair(OrbitalInfo(4, -4), 1));

        TwoM1.clear(); TwoM2.clear();
        TwoM2.push_back(1); TwoM2.push_back(1); TwoM2.push_back(1); TwoM2.push_back(1);

        Projection proj1(config1, TwoM1);
        Projection proj2(config2, TwoM2);
        *logstream << proj1.Name() << "* <- " << proj2.Name();
        many_body_operator.make_indirect_projection(proj1, indirects.left);
        many_body_operator.make_indirect_projection(proj2, indirects.right);

        int diffs = many_body_operator.GetProjectionDifferences<2>(indirects);
        *logstream << ((diffs > 0)? " = +< ": "= -< ") << indirects.left[0]->Name() << ", " << indirects.left[1]->Name()
                   << " |g| " << indirects.right[0]->Name() << ", " << indirects.right[1]->Name() << " >" << std::endl;
        EXPECT_EQ(2, diffs);
    }
}
