#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "Universal/Enums.h"
#include "Atom/MultirunOptions.h"
#include <iostream>
#include <set>

class Symmetry
{
    /** Storage for symmetry of Hamiltonian eigenstate (J, P), where
        J = total angular momentum and P = parity.
    */
public:
    Symmetry(int two_j, Parity parity):
        twoJ(two_j), P(parity)
    {}
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

/** Get set of requested symmetries from user input. */
std::set<Symmetry> ChooseSymmetries(const MultirunOptions& user_input);

#endif
