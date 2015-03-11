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
    EXPECT_DOUBLE_EQ(1./sqrt(2.), fabs(*it));
    it = ang4s5s.CSF_begin(1);
    EXPECT_DOUBLE_EQ(1./sqrt(2.), fabs(*it));
}

TEST(AngularDataTester, CountCSFs)
{
    // Create some RelativisticConfigurations and check the projections correspond to those expected.
    RelativisticConfiguration rconfig;

    //4d4
    rconfig.insert(std::make_pair(OrbitalInfo(4, -3), 2));
    rconfig.insert(std::make_pair(OrbitalInfo(4, 2), 2));

    int min_TwoJ = (rconfig.ParticleNumber()%2 == 0)? 0 : 1;
    AngularData ang(rconfig, min_TwoJ, min_TwoJ);

    int csf_count = 0;
    int proj_count_0 = ang.projection_size();

    for(int TwoJ = min_TwoJ; TwoJ <= rconfig.GetTwiceMaxProjection(); TwoJ += 2)
    {   AngularData ang(rconfig, TwoJ);
        csf_count += ang.GenerateCSFs(rconfig, TwoJ);
    }

    EXPECT_EQ(proj_count_0, csf_count);
}

TEST(AngularDataTester, HoleCSFs)
{
    // Create some electron only and electron-hole configurations and check that the
    // number of CSFs generated are equal.

    // i = Number of holes in 3d*
    for(int i = 0; i <= 2; i++)
    {
        RelativisticConfiguration rconfig_electron;
        RelativisticConfiguration rconfig_hole;

        // Electrons:
        // 3d8 4s1 = 3d(2+i) 3d*(6-i) 4s1
        if(6-i)
            rconfig_electron.insert(std::make_pair(OrbitalInfo(3, -3), 6-i));
        if(2+i)
            rconfig_electron.insert(std::make_pair(OrbitalInfo(3, 2), 2+i));
        rconfig_electron.insert(std::make_pair(OrbitalInfo(4, -1), 1));

        // Holes:
        // 3d-2 4s1 = 3d(-2+i) 3d*(-i) 4s1
        if(i)
            rconfig_hole.insert(std::make_pair(OrbitalInfo(3, -3), -i));
        if(i-2)
            rconfig_hole.insert(std::make_pair(OrbitalInfo(3, 2), i-2));
        rconfig_hole.insert(std::make_pair(OrbitalInfo(4, -1), 1));

        for(int TwoJ = 1; TwoJ <= 5; TwoJ += 2)
        {
            AngularData ang_electron(rconfig_electron, TwoJ);
            AngularData ang_hole(rconfig_hole, TwoJ);

            EXPECT_EQ(ang_electron.projection_size(), ang_hole.projection_size());

            int csf_count_electron = ang_electron.GenerateCSFs(rconfig_electron, TwoJ);
            int csf_count_hole = ang_hole.GenerateCSFs(rconfig_hole, TwoJ);

            EXPECT_EQ(csf_count_electron, csf_count_hole);
        }
    }
}

TEST(AngularDataTester, Iterators)
{
    // Create some RelativisticConfigurations and check the projections correspond to those expected.
    RelativisticConfiguration rconfig;

    //4d4
    rconfig.insert(std::make_pair(OrbitalInfo(4, -3), 3));
    rconfig.insert(std::make_pair(OrbitalInfo(4, 2), 2));

    int min_TwoJ = (rconfig.ParticleNumber()%2 == 0)? 0 : 1;
    pAngularDataLibrary library(new AngularDataLibrary(rconfig.ParticleNumber(), Symmetry(min_TwoJ, rconfig.GetParity()), min_TwoJ));

    rconfig.GetProjections(library);

    // Generate CSFs
    library->GenerateCSFs();

    // Norm check on CSFs
    ASSERT_LE(1, rconfig.NumCSFs());    // Make sure we're testing something

    std::vector<double> norm(rconfig.NumCSFs());
    auto it = rconfig.projection_begin(0);
    while(it != rconfig.projection_end(0))
    {
        auto csf_it = it.CSF_begin();
        int i = 0;
        while(csf_it != it.CSF_end())
        {   norm[i++] += (*csf_it) * (*csf_it);
            csf_it++;
        }

        it++;
    }
    for(int i = 0; i < norm.size(); i++)
        EXPECT_NEAR(1., norm[i], 1.e-6);
}

TEST(AngularDataTester, LadderOperator)
{
    // Create some RelativisticConfigurations and check the projections correspond to those expected.
    RelativisticConfiguration rconfig;

    //4s-1 4d4
    rconfig.insert(std::make_pair(OrbitalInfo(4, -1), -1));
    rconfig.insert(std::make_pair(OrbitalInfo(4, -3), 2));
    rconfig.insert(std::make_pair(OrbitalInfo(4, 2), 2));

    for(int parent_twoJ = 1; parent_twoJ <= 13; parent_twoJ += 2)
    {
        pAngularData parent(new AngularData(rconfig, parent_twoJ, parent_twoJ));
        int twoM = parent_twoJ - 2;

        while(twoM >= -parent_twoJ)
        {
            pAngularData child(new AngularData(rconfig, twoM));
            child->LadderLowering(rconfig, *parent);

            // Check normalisation
            const double* CSFs = child->GetCSFs();
            int num_csfs = child->NumCSFs();
            for(int csf_count = 0; csf_count < num_csfs; csf_count++)
            {
                double norm = 0.;
                for(int i = 0; i < child->projection_size(); i++)
                {
                    norm += CSFs[i * num_csfs + csf_count] * CSFs[i * num_csfs + csf_count];
                }

                EXPECT_NEAR(1.0, norm, 1.e-9);
            }

            parent = child;
            twoM = twoM - 2;
        }
    }

    // Verify AngularDataLibrary recursion
    int parent_twoJ = 7;
    pAngularData parent(new AngularData(rconfig, parent_twoJ, parent_twoJ));
    int twoM = parent_twoJ - 2;

    while(twoM >= 3)
    {
        pAngularData child(new AngularData(rconfig, twoM));
        child->LadderLowering(rconfig, *parent);
        
        parent = child;
        twoM = twoM - 2;
    }

    AngularDataLibrary lib(5, Symmetry(7, Parity::even), 3);
    pAngularData from_lib = lib[rconfig];
    lib.GenerateCSFs();

    // Compare CSFs
    ASSERT_EQ(from_lib->NumCSFs(), parent->NumCSFs());
    ASSERT_EQ(from_lib->projection_size(), parent->projection_size());
    const double* CSFs1 = from_lib->GetCSFs();
    const double* CSFs2 = parent->GetCSFs();
    for(int i = 0; i < from_lib->projection_size() * from_lib->NumCSFs(); i++)
        EXPECT_NEAR(CSFs1[i], CSFs2[i], 1.e-9);
}
