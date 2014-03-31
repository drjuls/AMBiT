#include "AngularData.h"
#include "RelativisticConfiguration.h"
#include "gtest/gtest.h"
#include "Include.h"
#include "HartreeFock/Core.h"

TEST(AngularDataTester, Projections)
{
    // Create some RelativisticConfigurations and check the projections correspond to those expected.
    RelativisticConfiguration rconfig;

    // Simple: 4s2
    rconfig.insert(std::make_pair(OrbitalInfo(4, -1), 2));

    // M = 0: expect [1, -1]
    AngularData ang4s2_0(rconfig, 0);
    ASSERT_EQ(1, ang4s2_0.projection_size());
    auto it = ang4s2_0.projection_begin();
    EXPECT_EQ(1, (*it)[0]);
    EXPECT_EQ(-1, (*it)[1]);

    // M = 1
    AngularData ang4s2_2(rconfig, 2);
    EXPECT_EQ(0, ang4s2_2.projection_size());

    // More complex: 4p 3d*2
    rconfig.clear();
    rconfig.insert(std::make_pair(OrbitalInfo(4, 1), 1));
    rconfig.insert(std::make_pair(OrbitalInfo(3, -3), 2));

    // TwoM = 3: First should be [-1, 3, 1]
    //           Last should be [1, 5, -3]
    AngularData ang3d24p1(rconfig, 3);
    ASSERT_EQ(4, ang3d24p1.projection_size());
    it = ang3d24p1.projection_begin();
    EXPECT_EQ(-1, (*it)[0]);
    EXPECT_EQ(3, (*it)[1]);
    EXPECT_EQ(1, (*it)[2]);

    // Point to last projection
    for(int count = 0; count < 3; count++)
        it++;
    EXPECT_EQ(1, (*it)[0]);
    EXPECT_EQ(5, (*it)[1]);
    EXPECT_EQ(-3, (*it)[2]);

    // Even more complex: 4f*3 6g*2
    rconfig.clear();
    rconfig.insert(std::make_pair(OrbitalInfo(6, -5), 2));
    rconfig.insert(std::make_pair(OrbitalInfo(4, -4), 3));

    // TwoM = 3: First should be [-1, -5, -7, 9, 7]
    //           Last should be [7, 5, 3, -3, -9]
    AngularData ang6g24f3_3(rconfig, 3);
    ASSERT_LE(1, ang6g24f3_3.projection_size());
    it = ang6g24f3_3.projection_begin();
    EXPECT_EQ(-1, (*it)[0]);
    EXPECT_EQ(-5, (*it)[1]);
    EXPECT_EQ(-7, (*it)[2]);
    EXPECT_EQ(9, (*it)[3]);
    EXPECT_EQ(7, (*it)[4]);

    // Point to last projection
    for(int count = 0; count < ang6g24f3_3.projection_size()-1; count++)
        it++;
    EXPECT_EQ(7, (*it)[0]);
    EXPECT_EQ(5, (*it)[1]);
    EXPECT_EQ(3, (*it)[2]);
    EXPECT_EQ(-3, (*it)[3]);
    EXPECT_EQ(-9, (*it)[4]);
}

TEST(AngularDataTester, CSFs)
{
    // Create some RelativisticConfigurations and check the projections correspond to those expected.
    RelativisticConfiguration rconfig;

    // Simple 4s2, M = 0: one CSF
    rconfig.insert(std::make_pair(OrbitalInfo(4, -1), 2));
    AngularData ang4s2_0(rconfig, 0, 0);
    ASSERT_EQ(1, ang4s2_0.NumCSFs());
    auto it = ang4s2_0.CSF_begin(0);
    EXPECT_DOUBLE_EQ(1., *it);

    // 4s5s M = 0, J = 0: two CSFs
    rconfig.clear();
    rconfig.insert(std::make_pair(OrbitalInfo(4, -1), 1));
    rconfig.insert(std::make_pair(OrbitalInfo(5, -1), 1));
    AngularData ang4s5s(rconfig, 0, 0);
    ASSERT_EQ(2, ang4s5s.projection_size());
    ASSERT_EQ(1, ang4s5s.NumCSFs());
    it = ang4s5s.CSF_begin(0);
    EXPECT_DOUBLE_EQ(1./sqrt(2.), fabs(*it++));
    EXPECT_DOUBLE_EQ(1./sqrt(2.), fabs(*it++));
}

TEST(AngularDataTester, CountCSFs)
{
    // Create some RelativisticConfigurations and check the projections correspond to those expected.
    RelativisticConfiguration rconfig;

    //4d4
    rconfig.insert(std::make_pair(OrbitalInfo(4, -3), 2));
    rconfig.insert(std::make_pair(OrbitalInfo(4, 2), 2));

    int min_TwoJ = (rconfig.ExcitationNumber()%2 == 0)? 0 : 1;
    AngularData ang(rconfig, min_TwoJ, min_TwoJ);

    int csf_count = 0;
    int proj_count_0 = ang.projection_size();

    for(int TwoJ = min_TwoJ; TwoJ <= rconfig.GetTwiceMaxProjection(); TwoJ += 2)
    {   AngularData ang(rconfig, TwoJ);
        csf_count += ang.GenerateCSFs(rconfig, TwoJ);
    }

    EXPECT_EQ(proj_count_0, csf_count);
}
