#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "Universal/Enums.h"
#include "Atom/MultirunOptions.h"
#include "HartreeFock/OrbitalInfo.h"
#include <iostream>
#include <set>

/** Symmetry provides a key to uniquely identify a particular Hamiltonian
    using an identifier for the Hamiltonian (J, parity).
    It is convertible to/from an int called Jpi, for which:
        Jpi >= 0:  TwoJ = Jpi,    Parity = even
        Jpi < 0:   TwoJ = -Jpi-1, Parity = odd
 */
class Symmetry
{
public:
    Symmetry(int Jpi): Jpi(Jpi) {}
    Symmetry(int two_j, Parity parity)
    {
        if(parity == Parity::even)
            Jpi = two_j;
        else
            Jpi = -two_j-1;
    }

    Symmetry(const std::string& twoJp)
    {
        std::stringstream ss(twoJp);
        ss >> Jpi;

        char c_parity;
        ss.get(c_parity);
        if(c_parity == 'o')
            Jpi = -Jpi-1;
    }

    Symmetry(const OrbitalInfo& orbital):
        Symmetry(orbital.TwoJ(), orbital.GetParity())
    {}

    Symmetry(const Symmetry& other): Jpi(other.Jpi) {}
    ~Symmetry() {}

    int GetTwoJ() const
    {
        if(Jpi < 0)
            return -Jpi-1;
        else
            return Jpi;
    }

    Parity GetParity() const
    {
        if(Jpi < 0)
            return Parity::odd;
        else
            return Parity::even;
    }

    double GetJ() const
    {   return double(GetTwoJ())/2.;
    }

    int GetJpi() const
    {   return Jpi;
    }

    /** Return string "<twoJ>.<P>" (e.g. "2.even"). */
    std::string GetString() const
    {
        std::string ret;
        std::stringstream angmom;
        angmom << GetTwoJ();
        ret = angmom.str();
        if(Jpi < 0)
            ret = ret + ".odd";
        else
            ret = ret + ".even";

        return ret;
    }

    /** Order first on parity P, then angular momentum J. */
    bool operator<(const Symmetry& other) const
    {
        if(Jpi >= 0)
        {
            if(other.Jpi < 0)
                return true;
            else
                return Jpi < other.Jpi;
        }
        else
        {
            if(other.Jpi >= 0)
                return false;
            else
                return Jpi > other.Jpi;
        }
    }

    bool operator==(const Symmetry& other) const
    {
        return Jpi == other.Jpi;
    }

    const Symmetry& operator=(const Symmetry& other)
    {
        Jpi = other.Jpi;
        return *this;
    }

protected:
    int Jpi;
};

#endif
