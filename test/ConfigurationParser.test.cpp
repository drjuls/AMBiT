#include "HartreeFock/ConfigurationParser.h"
#include "gtest/gtest.h"
#include "Include.h"
#include "Configuration/NonRelConfiguration.h"
#include "Configuration/RelativisticConfiguration.h"
#include "HartreeFock/Core.h"

using namespace Ambit;

TEST(ConfigurationParserTester, NonRelConfigs)
{
    std::string input("4d-1 5s2 5p4");

    NonRelConfiguration nrconfig(ConfigurationParser::ParseConfiguration<NonRelInfo, int>(input));
    EXPECT_EQ(-1, nrconfig.GetOccupancy(NonRelInfo(4, 2)));
    EXPECT_EQ(2, nrconfig.GetOccupancy(NonRelInfo(5, 0)));
    EXPECT_EQ(4, nrconfig.GetOccupancy(NonRelInfo(5, 1)));
    EXPECT_EQ(0, nrconfig.GetOccupancy(NonRelInfo(5, 2)));
}

TEST(ConfigurationParserTester, RelativisticConfigs)
{
    std::string input("4d+-1 5s2 5p2 5p+3");

    RelativisticConfiguration config(ConfigurationParser::ParseConfiguration<OrbitalInfo, int>(input));
    EXPECT_EQ(-1, config.GetOccupancy(OrbitalInfo(4, -3)));
    EXPECT_EQ(0, config.GetOccupancy(OrbitalInfo(4, 2)));
    EXPECT_EQ(2, config.GetOccupancy(OrbitalInfo(5, -1)));
    EXPECT_EQ(2, config.GetOccupancy(OrbitalInfo(5, 1)));
    EXPECT_EQ(3, config.GetOccupancy(OrbitalInfo(5, -2)));
}

TEST(ConfigurationParserTester, FractionalConfigs)
{
    std::string input("1s2 2s2 2p6 3s2 3p6 3d6");
    OccupationMap config;
    config = ConfigurationParser::ParseFractionalConfiguration(input);
    EXPECT_NEAR(2.0, config.GetOccupancy(OrbitalInfo(1, -1)), 1.e-6);
    EXPECT_NEAR(2.0, config.GetOccupancy(OrbitalInfo(2, 1)), 1.e-6);
    EXPECT_NEAR(4.0, config.GetOccupancy(OrbitalInfo(2, -2)), 1.e-6);
    EXPECT_NEAR(4, config.GetOccupancy(OrbitalInfo(3, 2)), 1.e-6);
    EXPECT_NEAR(2, config.GetOccupancy(OrbitalInfo(3, -3)), 1.e-6);

    input = "1s2 2s2 2p6 3s2 3p2";
    config = ConfigurationParser::ParseFractionalConfiguration(input);
    EXPECT_NEAR(2.0, config.GetOccupancy(OrbitalInfo(1, -1)), 1.e-6);
    EXPECT_NEAR(2.0, config.GetOccupancy(OrbitalInfo(2, 1)), 1.e-6);
    EXPECT_NEAR(4.0, config.GetOccupancy(OrbitalInfo(2, -2)), 1.e-6);
    EXPECT_NEAR(2.0, config.GetOccupancy(OrbitalInfo(3, 1)), 1.e-6);
    EXPECT_NEAR(0.0, config.GetOccupancy(OrbitalInfo(3, -2)), 1.e-6);

    input = "1s2 2s2 2p6 3s2 3p+2";
    config = ConfigurationParser::ParseFractionalConfiguration(input);
    EXPECT_NEAR(2.0, config.GetOccupancy(OrbitalInfo(1, -1)), 1.e-6);
    EXPECT_NEAR(2.0, config.GetOccupancy(OrbitalInfo(2, 1)), 1.e-6);
    EXPECT_NEAR(4.0, config.GetOccupancy(OrbitalInfo(2, -2)), 1.e-6);
    EXPECT_NEAR(0.0, config.GetOccupancy(OrbitalInfo(3, 1)), 1.e-6);
    EXPECT_NEAR(2.0, config.GetOccupancy(OrbitalInfo(3, -2)), 1.e-6);

    input = "1s2 2s2 2p6 3s2 3p+3.1 3p0.9";
    config = ConfigurationParser::ParseFractionalConfiguration(input);
    EXPECT_NEAR(2.0, config.GetOccupancy(OrbitalInfo(1, -1)), 1.e-6);
    EXPECT_NEAR(2.0, config.GetOccupancy(OrbitalInfo(2, 1)), 1.e-6);
    EXPECT_NEAR(4.0, config.GetOccupancy(OrbitalInfo(2, -2)), 1.e-6);
    EXPECT_NEAR(0.9, config.GetOccupancy(OrbitalInfo(3, 1)), 1.e-6);
    EXPECT_NEAR(3.1, config.GetOccupancy(OrbitalInfo(3, -2)), 1.e-6);
}

TEST(ConfigurationParserTester, BasisSize)
{
    std::string input("6s10pd4f");
    std::vector<int> basis(ConfigurationParser::ParseBasisSize(input));

    ASSERT_EQ(4, basis.size());
    EXPECT_EQ(6, basis[0]);
    EXPECT_EQ(10, basis[1]);
    EXPECT_EQ(10, basis[2]);
    EXPECT_EQ(4, basis[3]);

    input = " 6s 10pd4f ";
    basis = ConfigurationParser::ParseBasisSize(input);

    ASSERT_EQ(4, basis.size());
    EXPECT_EQ(6, basis[0]);
    EXPECT_EQ(10, basis[1]);
    EXPECT_EQ(10, basis[2]);
    EXPECT_EQ(4, basis[3]);
}

TEST(ConfigurationParserTester, Orbital)
{
    std::string input(" 5d");

    NonRelInfo nrinfo(ConfigurationParser::ParseOrbital<NonRelInfo>(input));
    EXPECT_TRUE(NonRelInfo(5, 2) == nrinfo);

    OrbitalInfo info(ConfigurationParser::ParseOrbital<OrbitalInfo>(input));
    EXPECT_TRUE(OrbitalInfo(5, 2) == info);
    EXPECT_FALSE(OrbitalInfo(5, -3) == info);

    input = "16p+";

    nrinfo = ConfigurationParser::ParseOrbital<NonRelInfo>(input);
    EXPECT_TRUE(NonRelInfo(16, 1) == nrinfo);

    info = ConfigurationParser::ParseOrbital<OrbitalInfo>(input);
    EXPECT_TRUE(OrbitalInfo(16, -2) == info);
    EXPECT_FALSE(OrbitalInfo(16, 1) == info);
}
