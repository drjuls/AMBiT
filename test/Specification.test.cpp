#include <array>
#include <string>

#include <gtest/gtest.h>
#include "Specification/Specification.h"

using namespace Ambit;

TEST(SpecificationTester, TransitionSectionAliases)
{
    for (const char* section : std::array{"Transitions", "Operators", "MatrixElements"}) {
        GlobalSpecification specification;
        const std::string input = std::string("[") + section + "/E1]\n"
            "MatrixElements = '1e:0->1e:1'\n"
            "RPA/BSpline/N = 42\n"
            "[" + section + "/HFS1]\n"
            "gOnI = 1.5\n"
            "[" + section + "/ResidualCoulomb]\n"
            "Scale = 2.0\n";

        ASSERT_TRUE(importSpecificationINIString(&specification, input).empty()) << section;
        EXPECT_TRUE(specification.transitions_do_e1) << section;
        ASSERT_TRUE(specification.transitions_e1.matrix_elements.has_value()) << section;
        EXPECT_EQ(specification.transitions_e1.matrix_elements->at(0), "1e:0->1e:1") << section;
        EXPECT_EQ(specification.transitions_e1.rpa_bspline_n, 42) << section;
        EXPECT_TRUE(specification.transitions_do_hfs1) << section;
        EXPECT_DOUBLE_EQ(specification.transitions_hfs1.g_on_i, 1.5) << section;
        EXPECT_TRUE(specification.transitions_do_residual_coulomb) << section;
        EXPECT_DOUBLE_EQ(specification.residual_coulomb.scale, 2.0) << section;
    }
}

TEST(SpecificationTester, TransitionSectionAliasesShareSettings)
{
    GlobalSpecification specification;
    const std::string input =
        "[Transitions/M1]\n"
        "Frequency = 0.25\n"
        "[Operators/M1]\n"
        "AllBelow = -1.0\n"
        "[MatrixElements/M1]\n"
        "Scale = 2.0\n";

    ASSERT_TRUE(importSpecificationINIString(&specification, input).empty());
    EXPECT_TRUE(specification.transitions_do_m1);
    EXPECT_DOUBLE_EQ(specification.transitions_m1.frequency, 0.25);
    ASSERT_TRUE(specification.transitions_m1.all_below.has_value());
    EXPECT_DOUBLE_EQ(*specification.transitions_m1.all_below, -1.0);
    EXPECT_DOUBLE_EQ(specification.transitions_m1.scale, 2.0);
}

TEST(SpecificationTester, TransitionSectionAliasesShareSettings)
{
    GlobalSpecification specification;
    const std::string input =
        "[Transitions]\n"
        "M1/Frequency = 0.25\n"
        "[Operators/M1]\n"
        "AllBelow = -1.0\n"
        "[MatrixElements]\n"
        "M1/Scale = 2.0\n";

    ASSERT_TRUE(importSpecificationINIString(&specification, input).empty());
    EXPECT_TRUE(specification.transitions_do_m1);
    EXPECT_DOUBLE_EQ(specification.transitions_m1.frequency, 0.25);
    ASSERT_TRUE(specification.transitions_m1.all_below.has_value());
    EXPECT_DOUBLE_EQ(*specification.transitions_m1.all_below, -1.0);
    EXPECT_DOUBLE_EQ(specification.transitions_m1.scale, 2.0);
}
