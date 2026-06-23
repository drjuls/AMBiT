#ifndef HAMILTONIAN_ID_H
#define HAMILTONIAN_ID_H

#include "Symmetry.h"
#include "RelativisticConfigList.h"
#include "Universal/Enums.h"
#include <boost/iterator/filter_iterator.hpp>

namespace Ambit
{
/** HamiltonianID provides a key to uniquely identify a particular Hamiltonian
    using an identifier for the Hamiltonian (e.g. J, parity).
 */
class HamiltonianID : public std::enable_shared_from_this<HamiltonianID>
{
public:
    HamiltonianID(int Jpi = 0): sym(Jpi) {}
    HamiltonianID(int two_J, Parity parity): sym(two_J, parity) {}
    HamiltonianID(const Symmetry& sym): sym(sym) {}
    HamiltonianID(const std::string& name);

    virtual bool operator<(const HamiltonianID& other) const
    {   return (sym < other.sym);
    }

    virtual bool operator>(const HamiltonianID& other) const
    {   return (other < *this);
    }

    virtual bool operator==(const HamiltonianID& other) const
    {   return (sym == other.sym);
    }

    double GetJ() const { return sym.GetJ(); }
    int GetTwoJ() const { return sym.GetTwoJ(); }
    Parity GetParity() const { return sym.GetParity(); }
    Symmetry GetSymmetry() const { return sym; }

    /** Unique short name for filenames, short printing.
        In this case, twoJ, followed by parity (e/o).
     */
    virtual std::string Name() const;

    /** Name for printing, in this case "J = <J>, P = <parity>". */
    virtual std::string Print() const;

    /** Read and write functions. */
    virtual void Write(FILE* fp) const;
    virtual void Read(FILE* fp);

    /** Clone object. Must be overridden by subclasses. */
    virtual std::shared_ptr<HamiltonianID> Clone() const
    {   return std::make_shared<HamiltonianID>(*this);
    }

protected:
    Symmetry sym;
};

typedef std::shared_ptr<HamiltonianID> pHamiltonianID;
typedef std::shared_ptr<const HamiltonianID> pHamiltonianIDConst;

inline std::ostream& operator<<(std::ostream& stream, const pHamiltonianID& hID)
{   return stream << hID->Print();
}

}
#endif
