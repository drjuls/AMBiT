#include "HamiltonianID.h"
#include "Include.h"
#include "Universal/MathConstant.h"
#include "NonRelConfiguration.h"

namespace Ambit
{
HamiltonianID::HamiltonianID(const std::string& name):
    sym(-1)
{
    int twoJ;
    std::stringstream ss(name);
    ss >> twoJ;

    Parity P;
    char c_parity;
    ss.get(c_parity);
    if(c_parity == 'o')
        P = Parity::odd;
    else
        P = Parity::even;

    sym = Symmetry(twoJ, P);
}

std::string HamiltonianID::Name() const
{
    std::string ret = itoa(sym.GetTwoJ());
    ret.append(LetterName(sym.GetParity()));

    return ret;
}

std::string HamiltonianID::Print() const
{
    std::string ret = "J = " + boost::lexical_cast<std::string>(sym.GetJ());
    ret += ", P = " + LowerName(sym.GetParity());

    return ret;
}

void HamiltonianID::Write(FILE* fp) const
{
    int two_j = sym.GetTwoJ();
    Parity P = sym.GetParity();
    file_err_handler->fwrite(&two_j, sizeof(int), 1, fp);
    file_err_handler->fwrite(&P, sizeof(Parity), 1, fp);
}

void HamiltonianID::Read(FILE* fp)
{
    int two_j;
    Parity P;
    file_err_handler->fread(&two_j, sizeof(int), 1, fp);
    file_err_handler->fread(&P, sizeof(Parity), 1, fp);
    sym = Symmetry(two_j, P);
}
}
