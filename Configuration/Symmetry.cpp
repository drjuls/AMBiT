#include "Symmetry.h"
#include "Include.h"

std::string Symmetry::GetString() const
{
    std::string ret;
    std::stringstream angmom;
    angmom << twoJ;
    ret = angmom.str();
    if(P == Parity::even)
        ret = ret + ".even";
    else
        ret = ret + ".odd";

    return ret;
}

bool Symmetry::operator<(const Symmetry& other) const
{
    if(P != other.P)
        return (P == Parity::even);
    else
        return (twoJ < other.twoJ);
}

bool Symmetry::operator==(const Symmetry& other) const
{
    if((P == other.P) && (twoJ == other.twoJ))
        return true;
    else
        return false;
}

const Symmetry& Symmetry::operator=(const Symmetry& other)
{
    P = other.P;
    twoJ = other.twoJ;
    
    return *this;
}
