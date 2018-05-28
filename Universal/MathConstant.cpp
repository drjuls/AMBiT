#include "Include.h"
#include "MathConstant.h"
#include <vector>
#include <algorithm>
#include <gsl/gsl_sf_coupling.h>

#ifdef AMBIT_USE_OPENMP
#include <omp.h>
#endif

MathConstant* MathConstant::Instance()
{
#ifdef AMBIT_USE_OPENMP
    // each OpenMP thread should get its own MathConstant Instance to maintain thread-safety when caching 3j/6j symbol values
    static std::unique_ptr<MathConstant[]> instances (new MathConstant[omp_get_max_threads()]());
    return &(instances[omp_get_thread_num()]);
#else
    // Obviously only return one instance if we're not using OpenMP
    static MathConstant instance;
    return &instance;
#endif
}

MathConstant::MathConstant():
    SpectroscopicNotation{"spdfghiklmnoqrtuv"}
{
    // Maximum largest stored twoJ for Wigner3J.
    // Set by consideration of key.
    MaxStoredTwoJ = 40;
    MSize = 2 * MaxStoredTwoJ + 1;

    // Initialise the densehash cache. The 3j key is strictly positive, so pick a negative key for the
    // empty bucket
    Symbols3j.set_empty_key(-1);
}

MathConstant::~MathConstant()
{}

double MathConstant::Pi() const { return 3.1415926535897932384626433832795; }

double MathConstant::BohrRadiusInFermi() const    { return 52917.7249; }
double MathConstant::HartreeEnergyIneV() const    { return 27.211383; }
double MathConstant::HartreeEnergyInInvCm() const { return 219474.631371; }
double MathConstant::SpeedOfLightSI() const       { return 299792458; }
double MathConstant::SpeedOfLightAU() const       { return 1./0.007297352932703; }
double MathConstant::InvCmInMHz() const           { return 29979.2458; }
double MathConstant::BohrRadiusSI() const         { return 5.2917720859e-11; }
double MathConstant::RydbergConstantSI() const    { return 10973731.6; }
double MathConstant::AtomicFrequencySI() const    { return 4.0 * Pi() * RydbergConstantSI() * SpeedOfLightSI(); }
double MathConstant::AtomicFrequencyMHz() const   { return RydbergConstantSI() * SpeedOfLightSI() * 2.e-6; }
double MathConstant::NuclearMagneton() const      { return 1.987131e-6; }
double MathConstant::Barn() const                 { return 1.e-28/(BohrRadiusSI() * BohrRadiusSI()); }
constexpr double MathConstant::ElectronMassInEV;

char MathConstant::GetSpectroscopicNotation(unsigned int l)  const
{
    if(l < SpectroscopicNotation.size())
        return SpectroscopicNotation[l];
    return 0;
}

int MathConstant::GetL(char spectroscopic_notation) const
{
    for(int L=0; L < SpectroscopicNotation.size(); L++)
    {
        if(SpectroscopicNotation[L] == spectroscopic_notation)
            return L;
    }

    return -1;
}

/** Clever method for finding wigner coefficients takes the log of the expression to convert
    the factorials into sums - that way the computer won't break. Also optimises by doing
    terms from the numerator and denominator in pairs so that we don't add and then subtract
    the same thing.
 */
double MathConstant::Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3) const
{
    /*  The whole expression is of the form
          ( a! b! ... c! )^(1/2)            ^(j1-j2-m3+k) n! m!
          ( ------------ )        x Sum (-1)             --------
          ( d! e! ... f! )           k                   p! q! r!
        Thus the expression has a part which is outside the summation ("outer") and
        a part inside the summation.
        Because of the summation we cannot take a logarithm of the whole thing. So we
        do each term in the sum separately, add them up, and multiply by the outer part.
     */
    if(m1 + m2 + m3 != 0.)
        return 0.;

    double test[9];
    test[0] = j1 + j2 - j3;
    test[1] = j3 + j1 - j2;
    test[2] = j2 + j3 - j1;
    test[3] = j1 + m1;
    test[4] = j1 - m1;
    test[5] = j2 + m2;
    test[6] = j2 - m2;
    test[7] = j3 + m3;
    test[8] = j3 - m3;

    // Check that all test[] are non-negative integers.
    unsigned int i;
    for(i=0; i<9; i++)
    {   if((test[i] < 0.) || (test[i] != floor(test[i])))
            return 0.;
    }

    double y[3];
    y[0] = j1 + j2 -j3;
    y[1] = j1 - m1;
    y[2] = j2 + m2;
    double z[3];
    z[0] = 0.;
    z[1] = j2 - j3 - m1;
    z[2] = j1 - j3 + m2;

    // Check that all y[] and z[] are integers.
    for(i=0; i<3; i++)
    {   if(y[i] != floor(y[i]))
            return 0.;

        if(z[i] != floor(z[i]))
            return 0.;
    }


    return gsl_sf_coupling_3j(int(2*j1), int(2*j2), int(2*j3), int(2*m1), int(2*m2), int(2*m3));
}

double MathConstant::Wigner6j(double j1, double j2, double j3, double j4, double j5, double j6) const
{
    /* The entire process is very similar to that used in calculating the 3j symbol. */
    if((j1 < 0.) || (j2 < 0.) || (j3 < 0.) || (j4 < 0.) || (j5 < 0.) || (j6 < 0.))
        return 0.;

    double test[4];
    test[0] = j1 + j2 + j3;
    test[1] = j1 + j5 + j6;
    test[2] = j4 + j2 + j6;
    test[3] = j4 + j5 + j3;

    double y[3];
    y[0] = j1 + j2 + j4 + j5;
    y[1] = j2 + j3 + j5 + j6;
    y[2] = j3 + j1 + j6 + j4;

    //unsigned int i;
    int i;

    for(i=0; i<4; i++)
    {   if(test[i] != floor(test[i]))
            return 0.;
    }
    for(i=0; i<3; i++)
    {   if(y[i] != floor(y[i]))
            return 0.;
    }
    
    double x[12];
    x[0] = j1 + j2 - j3;
    x[1] = j1 - j2 + j3;
    x[2] = - j1 + j2 + j3;
    x[3] = j1 + j5 - j6;
    x[4] = j1 - j5 + j6;
    x[5] = - j1 + j5 + j6;
    x[6] = j4 + j2 - j6;
    x[7] = j4 - j2 + j6;
    x[8] = - j4 + j2 + j6;
    x[9] = j4 + j5 - j3;
    x[10] = j4 - j5 + j3;
    x[11] = - j4 + j5 + j3;

    unsigned int j;
    for(j=0; j<12; j++)
    {   if((x[j] < 0.) || (x[j] != floor(x[j])))
            return 0.;
    }

    return gsl_sf_coupling_6j(int(2*j1), int(2*j2), int(2*j3), int(2*j4), int(2*j5), int(2*j6));
}

int MathConstant::HashWigner3j(int twoj1, int twoj2, int twoj3, int twom1, int twom2) const
{
    int key = twoj1 * MSize * MSize * 2 *(MaxStoredTwoJ + 1) * (MaxStoredTwoJ + 1)
            + twoj2 * MSize * MSize * 2 *(MaxStoredTwoJ + 1)
            + twoj3 * MSize * MSize     // twoj3 can be as large as 2 * MaxStoredTwoJ (when k = j1 + j2)
            + twom1 * MSize             // m1 >= 0 && m1 > m2
            + twom2 + MaxStoredTwoJ;
    return key;
}

double MathConstant::Electron3j(int twoj1, int twoj2, int k, int twom1, int twom2)
{
    // Check for physical correctness:
    //  - triangle condition (j1, j2, k)
    //  - k >= abs(q)
    int twok = 2 * k;
    int twoq = - twom1 - twom2;
    if((twok > twoj1 + twoj2) || (abs(twoj1 - twoj2) > twok)
       || (twok < abs(twoq)) || (twoj1 < abs(twom1)) || (twoj2 < abs(twom2)))
        return 0.;

    // Sort such that j1 >= j2, m1 >=0, and if (j1 == j2) then m1 >= m2.
    // Keep track of sign changes using boolean (false -> take negative)
    bool sign_swap = ((twoj1 + twoj2)/2 + k)%2 == 1;
    bool sign = true;

    // Sort so that (j1 >= j2) for lookup
    if(twoj1 < twoj2)
    {
        std::swap(twoj1, twoj2);
        std::swap(twom1, twom2);
        if(sign_swap) sign = !sign;
    }

    // (m1 >= 0)
    if(twom1 < 0)
    {   twom1 = -twom1; twom2 = -twom2;
        twoq = -twoq;
        if(sign_swap) sign = !sign;
    }

    // if(j1 == j2), then (m1 >= m2)
    if((twoj1 == twoj2) && (twom1 < twom2))
    {
        std::swap(twom1, twom2);
        if(sign_swap) sign = !sign;
    }

    if(twoj1 > MaxStoredTwoJ)
    {   
        if(sign)
            return gsl_sf_coupling_3j(twoj1, twoj2, 2*k, twom1, twom2, twoq);
        else
            return -gsl_sf_coupling_3j(twoj1, twoj2, 2*k, twom1, twom2, twoq);
    }

    size_t key = HashWigner3j(twoj1, twoj2, twok, twom1, twom2);
    google::dense_hash_map<int, double>::const_iterator it = Symbols3j.find(key);

    if(it != Symbols3j.end())
    {   if(sign)
            return it->second;
        else
            return -it->second;
    }
    else
    {   
        double value = gsl_sf_coupling_3j(twoj1, twoj2, 2*k, twom1, twom2, twoq);
        Symbols3j[key] = value;
        
        if(sign)
            return value;
        else
            return -value;
    }
}

double MathConstant::Electron3j(int twoj1, int twoj2, int k)
{
    // In this case
    //  ( j1   j2   k ) = (-1)^(j1 + j2 + k) (  j2    j1   k ) = ( j2   j1   k )
    //  ( 1/2 -1/2  0 )                      ( -1/2  1/2   0 )   ( 1/2 -1/2  0 )
    // So we don't need to worry about sign.
    // Therefore we can cut out some tests by asserting j1 >= j2

    if(twoj1 >= twoj2)
        return Electron3j(twoj1, twoj2, k, 1, -1);
    else
        return Electron3j(twoj2, twoj1, k, 1, -1);
}

double MathConstant::Wigner3j(int twoj1, int twoj2, int twoj3, int twom1, int twom2)
{
    int twom3 = -twom1 -twom2;

    // Check if it has the form of Electron3j (two half-integer, one integer)
    if(twoj1%2 && twoj2%2 && !twoj3%2)
        return Electron3j(twoj1, twoj2, twoj3/2, twom1, twom2);
    else if(twoj1%2 && !twoj2%2 && twoj3%2)
        return Electron3j(twoj3, twoj1, twoj2/2, twom3, twom1);
    else if(!twoj1%2 && twoj2%2 && twoj3%2)
        return Electron3j(twoj2, twoj3, twoj1/2, twom2, twom3);

    // Otherwise must be all integer
    // Check for physical correctness:
    //  - triangle condition (j1, j2, k)
    //  - k >= abs(q)
    if((twoj3 > twoj1 + twoj2) || (abs(twoj1 - twoj2) > twoj3)
       || (twoj1 < abs(twom1)) || (twoj2 < abs(twom2)) || (twoj3 < abs(twom3)))
        return 0.;

    // Arrange so that j1 >= j2 >= j3, m1 >= 0 and m1 > m2 > m3 where possible
    // Keep track of sign changes using boolean (false -> take negative)
    bool sign_swap = ((twoj1 + twoj2 + twoj2)/2)%2 == 1;
    bool sign = true;

    // Make sure j1 is the largest
    if(twoj2 > twoj1 && twoj2 >= twoj3)
    {
        std::swap(twoj1, twoj2);
        std::swap(twom1, twom2);
        if(sign_swap) sign = !sign;
    }
    else if(twoj3 > twoj1 && twoj3 > twoj2)
    {
        std::swap(twoj1, twoj3);
        std::swap(twom1, twom3);
        if(sign_swap) sign = !sign;
    }

    // Sort so that (j2 >= j3)
    if(twoj2 < twoj3)
    {
        std::swap(twoj2, twoj3);
        std::swap(twom2, twom3);
        if(sign_swap) sign = !sign;
    }

    if(twoj1 > MaxStoredTwoJ)
    {
        if(sign)
            return gsl_sf_coupling_3j(twoj1, twoj2, twoj3, twom1, twom2, twom3);
        else
            return -gsl_sf_coupling_3j(twoj1, twoj2, twoj3, twom1, twom2, twom3);
    }

    // (m1 >= 0)
    if(twom1 < 0)
    {   twom1 = -twom1; twom2 = -twom2; twom3 = -twom3;
        if(sign_swap) sign = !sign;
    }

    // if(j1 == j2 == j3), then (m1 >= m2  && m1 >= m3)
    if(twoj1 == twoj3)
    {
        if(twom2 > twom1 && twom2 >= twom3)
        {
            std::swap(twom1, twom2);
            if(sign_swap) sign = !sign;
        }
        else if(twom3 > twom1 && twom3 > twom2)
        {
            std::swap(twom1, twom3);
            if(sign_swap) sign = !sign;
        }
    }
    // if(j1 == j2), then (m1 >= m2)
    else if((twoj1 == twoj2) && (twom1 < twom2))
    {
        std::swap(twom1, twom2);
        if(sign_swap) sign = !sign;
    }

    // if(j2 == j3), then (m2 >= m3)
    if((twoj2 == twoj3) && twom2 < twom3)
    {
        std::swap(twom2, twom3);
        if(sign_swap) sign = !sign;
    }

    size_t key = HashWigner3j(twoj1, twoj2, twoj3, twom1, twom2);
    google::dense_hash_map<int, double>::const_iterator it = Symbols3j.find(key);

    if(it != Symbols3j.end())
    {   if(sign)
            return it->second;
        else
            return -it->second;
    }
    else
    {   double value = gsl_sf_coupling_3j(twoj1, twoj2, twoj3, twom1, twom2, twom3);
        Symbols3j[key] = value;
        
        if(sign)
            return value;
        else
            return -value;
    }
}

unsigned int MathConstant::GetStorageSize() const
{
    return Symbols3j.size();
}

void MathConstant::Reset()
{
    Symbols3j.clear();
}

double MathConstant::SphericalTensorReducedMatrixElement(int kappa1, int kappa2, int rank)
{
    int l1 = 0;
    int l2 = 0;
    int two_j1 = 0;
    int two_j2 = 0;

    if(kappa1 < 0)
    {
        two_j1 = - 2 * kappa1 - 1;
        l1 = -kappa1 - 1;
    }
    else
    {
        two_j1 = 2 * kappa1 - 1;
        l1 = kappa1;
    }

    if(kappa2 < 0)
    {
        two_j2 = - 2 * kappa2 - 1;
        l2 = -kappa2 - 1;
    }
    else
    {
        two_j2 = 2 * kappa2 - 1;
        l2 = kappa2;
    }

    if((l1 + rank + l2)%2 != 0)
    {
        return 0.0;
    }

    double value = MathConstant::Instance()->Electron3j(two_j1, two_j2, rank);
    if(value)
    {
        value *= minus_one_to_the_power(two_j1 + (two_j2+1)/2 + rank);
        value *= sqrt((two_j1 + 1.0) * (two_j2 + 1.0));
    }

    return value;
}

unsigned int MathConstant::nChoosek(unsigned int n, unsigned int k) const
{
    if(k > n)
        return 0;
    if(k * 2 > n)
        k = n-k;
    if(k == 0)
        return 1;

    unsigned int result = n;
    for(int i = 2; i <= k; ++i )
    {   result *= (n-i+1);
        result /= i;
    }
    return result;
}
