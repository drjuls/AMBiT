#include "Include.h"
#include "Constant.h"
#include <vector>
#include <algorithm>

const double Constant::Pi = 3.1415926535897932384626433832795;
double Constant::Alpha = 7.297352932703e-3;
double Constant::AlphaSquared = Constant::Alpha * Constant::Alpha;
const double Constant::NucleonElectronMassRatio = 1822.888;

const double Constant::AtomicToFermi = 52917.7249;
const double Constant::HartreeEnergy_eV = 27.211383;
const double Constant::HartreeEnergy_cm = 219474.631371;
const double Constant::SpeedOfLight = 299792458;
const double Constant::InvCmToMHz = 29979.2458;

const char Constant::SpectroscopicNotation[10] 
    = {'s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l'};

/** Clever method for finding wigner coefficients takes the log of the expression to convert
    the factorials into sums - that way the computer won't break. Also optimises by doing
    terms from the numerator and denominator in pairs so that we don't add and then subtract
    the same thing.
 */
double LogFactorialFraction(unsigned int num, unsigned int denom);
double Constant::Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3)
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

    std::vector<unsigned int> outer_num(9), outer_den(9);   // Numerator terms and denominator
                                                            // terms of outer part.
    // Check that all test[] are non-negative integers. Save the integers.
    unsigned int i;
    for(i=0; i<9; i++)
    {   if((test[i] < 0.) || (test[i] != floor(test[i])))
            return 0.;
        outer_num[i] = (unsigned int)floor(test[i]);
    }

    double y[3];
    y[0] = j1 + j2 -j3;
    y[1] = j1 - m1;
    y[2] = j2 + m2;
    double z[3];
    z[0] = 0.;
    z[1] = j2 - j3 - m1;
    z[2] = j1 - j3 + m2;

    // Check that all y[] and z[] are integers. Save the integers.
    std::vector<int> ma(3), na(3);
    for(i=0; i<3; i++)
    {   if(y[i] != floor(y[i]))
            return 0.;
        ma[i] = (int)(y[i]);

        if(z[i] != floor(z[i]))
            return 0.;
        na[i] = (int)(z[i]);
    }

    std::sort(ma.begin(), ma.end());
    std::sort(na.begin(), na.end());
    // All elements of ma should be non-negative
    if(ma[0] < 0)
        return 0.;
    // Each element of ma should be larger than all elements of na.
    if(ma[0] < na[2])
        return 0.;

    // Calculate denominator terms in outer part
    outer_den[0] = (unsigned int)(j1 + j2 + j3 + 1.);
    unsigned int j;
    for(j=1; j<=2; j++)
    {   outer_den[j] = outer_den[j+2] = ma[j] - ma[0];
        outer_den[j+4] = outer_den[j+6] = na[2] - na[j-1];
    }

    // Calculate outer part
    std::sort(outer_den.begin(), outer_den.end());
    std::sort(outer_num.begin(), outer_num.end());

    double outer = 0.;
    for(j=0; j<9; j++)
    {   outer += LogFactorialFraction(outer_num[j], outer_den[j]);
    }
    outer = outer/2;    // Outer part is to the power 1/2.

    // Calculate summation part
    double total = 0.;
    for(int k = na[2]; k <= ma[0]; k++)
    {   double term = 0.;
        term += LogFactorialFraction(1, ma[0] - k);
        term += LogFactorialFraction(ma[1] - ma[0], ma[1] - k);
        term += LogFactorialFraction(ma[2] - ma[0], ma[2] - k);
        term += LogFactorialFraction(0, k - na[2]);
        term += LogFactorialFraction(na[2] - na[0], k - na[0]);
        term += LogFactorialFraction(na[2] - na[1], k - na[1]);
        double L = k + j1 -j2 - m3;
        total = total + pow(-1., L) * exp(term);
    }
    total = exp(outer) * total;
    
    return total;
}

double Constant::Wigner6j(double j1, double j2, double j3, double j4, double j5, double j6)
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

    std::vector<unsigned int> ma(4), na(3);
    unsigned int i;
    for(i=0; i<4; i++)
    {   if(test[i] != floor(test[i]))
            return 0.;
        ma[i] = (unsigned int)test[i];
    }
    for(i=0; i<3; i++)
    {   if(y[i] != floor(y[i]))
            return 0.;
        na[i] = (unsigned int)y[i];
    }

    std::sort(ma.begin(), ma.end());
    std::sort(na.begin(), na.end());

    if(ma[3] > na[1])
        return 0.;

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

    std::vector<unsigned int> outer_den(13), outer_num(13);
    unsigned int j;
    for(j=0; j<12; j++)
    {   if((x[j] < 0.) || (x[j] != floor(x[j])))
            return 0.;
        outer_num[j] = (unsigned int)x[j];
    }
    outer_num[12] = ma[3] + 1;
    for(j=0; j<3; j++)
    {   outer_den[j] = ma[j] + 1;
        outer_den[j+3] = outer_den[j+6] = ma[3] - ma[j];
    }
    for(j=0; j<2; j++)
    {   outer_den[j+9] = outer_den[j+11] = na[j+1] - na[0];
    }
    std::sort(outer_den.begin(), outer_den.end());
    std::sort(outer_num.begin(), outer_num.end());

    double outer = 0.;
    for(j=0; j<13; j++)
    {   outer += LogFactorialFraction(outer_num[j], outer_den[j]);
    }
    outer = outer/2;

    double total = 0.;
    for(unsigned int k=ma[3]; k<=na[0]; k++)
    {   double term = 0.;
        term += LogFactorialFraction(k+1, ma[3]+1);
        term += LogFactorialFraction(ma[3] - ma[0], k - ma[0]);
        term += LogFactorialFraction(ma[3] - ma[1], k - ma[1]);
        term += LogFactorialFraction(ma[3] - ma[2], k - ma[2]);
        term += LogFactorialFraction(0, k - ma[3]);
        term += LogFactorialFraction(0, na[0] - k);
        term += LogFactorialFraction(na[1] - na[0], na[1] - k);
        term += LogFactorialFraction(na[2] - na[0], na[2] - k);
        total = total + pow(-1., double(k)) * exp(term);
    }
    total = exp(outer) * total;

    return total;
}

/** Calculate the logarithm of a fraction where the numerator and denominator are factorials
     log( n!/d! )
 */
double LogFactorialFraction(unsigned int num, unsigned int denom)
{
    double total = 0.;
    if(num > denom)
    {   for(unsigned int i=denom+1; i<=num; i++)
            total = total + log(double(i));
    }
    else if(denom > num)
    {   for(unsigned int i=num+1; i<=denom; i++)
            total = total - log(double(i));
    }
    return total;
}

const unsigned int Constant::MaxStoredTwoJ = 11;    // Up to h-shell electrons
std::map<int, double> Constant::Symbols3j;

double Constant::Electron3j(unsigned int twoj1, unsigned int twoj2, unsigned int k, int twom1, int twom2)
{
    double sign = 1.;

    // Sort so that (j1 >= j2) for lookup
    if(twoj1 < twoj2)
    {   unsigned int tempj;
        tempj = twoj1; twoj1 = twoj2; twoj2 = tempj;
        int tempm;
        tempm = twom1; twom1 = twom2; twom2 = tempm;
        if(((twoj1 + twoj2)/2 + k)%2 == 1)
            sign = -sign;
    }
    else if((twoj1 == twoj2) && (twom1 < twom2))
    {   int tempm;
        tempm = twom1; twom1 = twom2; twom2 = tempm;
        if((1 + k)%2 == 1)
            sign = -sign;
    }

    if(twoj1 > MaxStoredTwoJ)
    {   double q = double(- twom1 - twom2)/2.;
        return sign * Wigner3j(double(twoj1)/2., double(twoj2)/2., double(k), double(twom1)/2., double(twom2)/2., q);
    }
    else if((2*k > twoj1 + twoj2) || (abs(int(twoj1) - int(twoj2)) > 2*int(k)))
    {   return 0.;
    }

    const int MSize = 2 * MaxStoredTwoJ + 1;

    int key = int(twoj1) * MSize * MSize * (MaxStoredTwoJ + 1) * MaxStoredTwoJ
            + int(twoj2) * MSize * MSize * (MaxStoredTwoJ + 1)
            + int(k)     * MSize * MSize
                + twom1  * MSize
                + twom2;
    std::map<int, double>::const_iterator it = Symbols3j.find(key);

    if(it != Symbols3j.end())
    {   return sign * it->second;
    }
    else
    {   double q = double(- twom1 - twom2)/2.;
        double value = Wigner3j(double(twoj1)/2., double(twoj2)/2., double(k), double(twom1)/2., double(twom2)/2., q);
        Symbols3j[key] = value;
        return sign * value;
    }
}