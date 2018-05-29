#ifndef LEVEL_H
#define LEVEL_H

#include "Symmetry.h"
#include "RelativisticConfigList.h"
#include "Universal/Enums.h"
#include <boost/iterator/filter_iterator.hpp>

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

/** A Level is an eigenstate of the Hamiltonian with eigenvector of length NumCSFs.
    The symmetry of the eigenstate is given by the stored HamiltonianID.
 */
class Level : public std::enable_shared_from_this<Level>
{
public:
    /** Initialise Level with energy eigenvalue, eigenvector, and pointer to relativistic configurations (CSFs).
        PRE: length(csf_eigenvector) == configlist->NumCSFs() == numCSFs (if supplied).
     */
    Level(const double& energy, const std::vector<double>& csf_eigenvector, pHamiltonianID hamiltonian_id, const double& gFactor):
        eigenvalue(energy), eigenvector(csf_eigenvector), hamiltonian(hamiltonian_id), gFactor(gFactor) {}
    Level(const double& energy, const double* csf_eigenvector, pHamiltonianID hamiltonian_id, unsigned int numCSFs);

    double GetEnergy() const { return eigenvalue; }
    void SetEnergy(double energy) { eigenvalue = energy; }

    double GetgFactor() const { return gFactor; }
    void SetgFactor(double g_factor) { gFactor = g_factor; }

    const std::vector<double>& GetEigenvector() const { return eigenvector; }
    unsigned int GetEigenvectorLength() const { return eigenvector.size(); }     //!< Eigenvector length = NumCSFs

    Parity GetParity() const { return hamiltonian->GetParity(); }
    unsigned int GetTwoJ() const { return hamiltonian->GetTwoJ(); }

    pHamiltonianID GetHamiltonianID() { return hamiltonian; }
    pHamiltonianIDConst GetHamiltonianID() const { return hamiltonian; }

protected:
    double eigenvalue;
    std::vector<double> eigenvector;        // length of eigenvector = hamiltonian->configs->NumCSFs()
    pHamiltonianID hamiltonian = nullptr;   // identifier for the Hamiltonian that *this is an eigenvalue of
    double gFactor = std::numeric_limits<double>::quiet_NaN();
};

typedef std::shared_ptr<Level> pLevel;
typedef std::shared_ptr<const Level> pLevelConst;

#endif
