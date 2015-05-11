#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "Universal/Enums.h"
#include "Atom/MultirunOptions.h"
#include <iostream>
#include <set>

/** Symmetry provides a key to uniquely identify a particular Hamiltonian
    using an identifier for the Hamiltonian (J, parity).
 */
class Symmetry
{
public:
    Symmetry(int kappa)
    {   twoJ = 2*abs(kappa)-1;
        int L = (kappa > 0)? kappa : -kappa-1;
        P = (L%2 == 0)? Parity::even : Parity::odd;
    }
    Symmetry(int two_j, Parity parity):
        twoJ(two_j), P(parity)
    {}
    Symmetry(const std::string& twoJp);
    Symmetry(const Symmetry& other):
        twoJ(other.twoJ), P(other.P)
    {}
    ~Symmetry() {}

    int GetTwoJ() const
    {   return twoJ;
    }

    Parity GetParity() const
    {   return P;
    }

    double GetJ() const
    {   return double(twoJ)/2.;
    }

    int GetKappa() const
    {
        int kappa = (twoJ + 1)/2;
        if(kappa%2 == 1)
            kappa = -kappa;
        if(P == Parity::odd)
            kappa = -kappa;

        return kappa;
    }

    /** Return string "<twoJ>.<P>" (e.g. "2.even"). */
    std::string GetString() const;

    /** Order first on parity P, then angular momentum J. */
    bool operator<(const Symmetry& other) const;
    bool operator==(const Symmetry& other) const;
    const Symmetry& operator=(const Symmetry& other);

protected:
    int twoJ;
    Parity P;
};

#endif
