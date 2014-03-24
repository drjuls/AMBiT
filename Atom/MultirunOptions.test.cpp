#include "Atom/MultirunOptions.h"
#include "Atom/GetPot"
#include "gtest/gtest.h"
#include "Include.h"

TEST(MultirunOptionsTester, MultirunTest)
{
    std::string inputFileName("template.input");
    MultirunOptions opt(inputFileName.c_str(), "//", "\n", ",");

    // Check AlphaSquaredVariation (the multirun double)
    EXPECT_EQ(3, opt.GetNumRuns());
    opt.SetRun(2);
    EXPECT_DOUBLE_EQ(0.1, opt("AlphaSquaredVariation", 0.0));

    // Check non-multirun doubles act normally
    EXPECT_DOUBLE_EQ(2.3, opt("NuclearThickness", 0.0));
    EXPECT_DOUBLE_EQ(0.000001, opt("Lattice/StartPoint", 0.0));

    EXPECT_STREQ("test", opt("ID", "").c_str());
}

TEST(MultirunOptionsTester, AbsorbTest)
{
    MultirunOptions opt;

    std::string inputFileName("template.input");
    MultirunOptions opt2(inputFileName.c_str(), "//", "\n", ",");

    opt.absorb(opt2);

    // Check AlphaSquaredVariation (the multirun double)
    EXPECT_EQ(3, opt.GetNumRuns());
    opt.SetRun(0);
    EXPECT_DOUBLE_EQ(-0.1, opt("AlphaSquaredVariation", 0.0));

    // Check non-multirun doubles act normally
    EXPECT_DOUBLE_EQ(2.3, opt("NuclearThickness", 0.0));
//    EXPECT_DOUBLE_EQ(0.000001, opt("Lattice/StartPoint", 0.0));
//    EXPECT_EQ(true, opt.search("Basis/--bspline-basis"));
}

TEST(MultirunOptionsTester, SearchTest)
{
    std::string inputFileName("template.input");
    MultirunOptions opt(inputFileName.c_str(), "//", "\n", ",");

    // Check search works
    EXPECT_EQ(true, opt.search("Basis/--bspline-basis"));
    EXPECT_EQ(false, opt.search("--simply-false"));
}

TEST(MultirunOptionsTester, VectorTest)
{
    std::string inputFileName("template.input");
    MultirunOptions opt(inputFileName.c_str(), "//", "\n", ",");

    ASSERT_EQ(2, opt.vector_variable_size("CI/LeadingConfigurations"));
    EXPECT_STREQ("4s2", opt("CI/LeadingConfigurations", "", 0).c_str());
}

//TEST(MultirunOptionsTester, UFOTest)
//{
//    std::string inputFileName("template.input");
//    MultirunOptions opt(inputFileName.c_str(), "//", "\n", ",");
//
//    EXPECT_EQ(18, opt("HF/N", 0));
//
//    opt.set_prefix("HF/");
//    EXPECT_EQ(18, opt("N", 1));
//
//    std::vector<std::string> ufo = opt.unidentified_variables();
//    std::vector<std::string>::const_iterator it = ufo.begin();
//    while(it != ufo.end())
//    {   *outstream << *it << std::endl;
//        it++;
//    }
//}
