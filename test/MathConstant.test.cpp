#include "Universal/MathConstant.h"
#include "Universal/FornbergDifferentiator.h"
#include "Include.h"
#include "gtest/gtest.h"

using namespace Ambit;

TEST(ConstantTester, Electron3j)
{
    MathConstant* constant = MathConstant::Instance();

    constant->Reset();

    // Sanity check
    EXPECT_DOUBLE_EQ(-1.0/std::sqrt(10.0), constant->Electron3j(5, 3, 1, 3, -1));
    EXPECT_DOUBLE_EQ(1.0/std::sqrt(10.0), constant->Electron3j(3, 5, 1, -1, 3));

    // Should only have used one storage unit
    EXPECT_EQ(1, constant->GetStorageSize());

    // Fill 3j storage and check number of elements
    int twoj1, twoj2, k;
    int twom1, twom2;

    double electron3j_value, wigner3j_value;

    // Max j is 10 for this test.
    unsigned int maxtwoj = 20;
    for(twoj1 = 1; twoj1 <= maxtwoj; twoj1+=2)
    {
        // j1 and j2 are both integer or half-integer, hence the mod 2
        for(twoj2 = twoj1%2; twoj2 <= maxtwoj; twoj2+=2)
        {
            for(k = 0; k <= 2*maxtwoj; k++)
            {
                // Span all m values - let MathConstant::Electron3j() sort it out
                for(twom1 = 1; twom1 <= twoj1; twom1+=2)
                {
                    for(twom2 = -twoj2; twom2 <= twoj2; twom2+=2)
                    {
                        electron3j_value = constant->Electron3j(twoj1, twoj2, k, twom1, twom2);
                        wigner3j_value = constant->Wigner3j(double(twoj1)/2., double(twoj2)/2., double(k),
                                                            double(twom1)/2., double(twom2)/2., double(-twom1-twom2)/2.);

                        if(std::fabs(wigner3j_value) > 1.e-12)
                            EXPECT_NEAR(wigner3j_value, electron3j_value, std::fabs(1.e-10*electron3j_value));
                    }
                }
            }
        }
    }

    // Size of storage here calculated using Mathematica
    EXPECT_EQ(28248, constant->GetStorageSize());

    // Finally, check three-parameter version
    // ( 3/2  7/2 4 ) = sqrt(5/7)/6
    // ( 1/2 -1/2 0 )
    EXPECT_DOUBLE_EQ(std::sqrt(5./7.)/6., constant->Electron3j(3, 7, 4));

    constant->Reset();
    EXPECT_EQ(0, constant->GetStorageSize());
}

TEST(ConstantTester, SpectroscopicNotation)
{
    MathConstant* constant = MathConstant::Instance();
    EXPECT_EQ('d', constant->GetSpectroscopicNotation(2));
    EXPECT_EQ(3, constant->GetL('f'));
}

TEST(ConstantTester, Wigner6j)
{
    EXPECT_DOUBLE_EQ(-0.1, MathConstant::Instance()->Wigner6j(2.5, 1.5, 2., 1.5, 1.5, 2.));
    EXPECT_DOUBLE_EQ(-std::sqrt(2./7.)/5., MathConstant::Instance()->Wigner6j(2.5, 1.5, 2., 1.5, 1.5, 3.));
    EXPECT_DOUBLE_EQ(0.0, MathConstant::Instance()->Wigner6j(2.5, 1.5, 2., 1.5, 2.5, 4.));
}

TEST(FornbergTester, Sine)
{
    pLattice lattice = std::make_shared<Lattice>(1000, 1.e-6, 50);
    FornbergDifferentiator diff(lattice, 7, true);
    const double* R = lattice->R();

    std::vector<double> sine(1000);
    std::vector<double> cos(1000);
    for(int i = 0; i < sine.size(); i++)
    {
        sine[i] = std::sin(R[i]);
        cos[i] = std::cos(R[i]);
    }

    // First derivative
    std::vector<double> dsine(1000);
    diff.GetDerivative(sine, dsine);
    for(int i = 0; i < sine.size(); i++)
    {
        EXPECT_NEAR(cos[i], dsine[i], 1.e-6);
    }

    std::vector<double> d2sine(1000);
    diff.GetSecondDerivative(sine, d2sine);

    for(int i = 0; i < sine.size(); i++)
    {
        EXPECT_NEAR(-sine[i], d2sine[i], 2.e-5);
    }
}
