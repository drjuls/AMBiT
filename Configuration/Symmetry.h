#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <iostream>

#ifndef PARITY_ENUM
#define PARITY_ENUM
    enum Parity { even, odd };
#endif

class Symmetry
{
    /** Storage for symmetry of Hamiltonian eigenstate (J, P), where
        J = total angular momentum and P = parity.
    */
public:
    Symmetry(unsigned int two_j, Parity parity):
        twoJ(two_j), P(parity)
    {}
    Symmetry(const Symmetry& other):
        twoJ(other.twoJ), P(other.P)
    {}
    ~Symmetry() {}

    unsigned int GetTwoJ() const
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
    unsigned int twoJ;
    Parity P;
};

#endif
