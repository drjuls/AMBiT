#include "Atom/MultirunOptions.h"
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
    //EXPECT_DOUBLE_EQ(0.000001, opt("Lattice/StartPoint", 0.0));
    //EXPECT_EQ(true, opt.search("Basis/--bspline-basis"));
}
