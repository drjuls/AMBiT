#ifndef HAMILTONIAN_TYPES_H
#define HAMILTONIAN_TYPES_H

#include "Configuration/Level.h"

/** SingleOrbitalID defines a Hamiltonian of size 1:
    just the one configuration from one orbital.
    They are distinguishable by PQN().
 */
class SingleOrbitalID : public HamiltonianID
{
public:
    SingleOrbitalID(): HamiltonianID(), pqn(0) {}
    SingleOrbitalID(const OrbitalInfo& info): HamiltonianID(Symmetry(info)), pqn(info.PQN()) {}
    SingleOrbitalID(const SingleOrbitalID& other): HamiltonianID(other), pqn(other.pqn) {}

    virtual bool operator<(const HamiltonianID& other) const override
    {
        if(sym < other.GetSymmetry())
            return true;
        else if (other.GetSymmetry() < sym)
            return false;
        else
        {   const SingleOrbitalID* other_soid = dynamic_cast<const SingleOrbitalID*>(&other);
            if(other_soid)
                return pqn < other_soid->pqn;
            else
                return false;
        }
    }
    
    virtual bool operator==(const HamiltonianID& other) const
    {
        const SingleOrbitalID* other_soid = dynamic_cast<const SingleOrbitalID*>(&other);
        if(other_soid)
            return (sym == other_soid->sym) && (pqn == other_soid->pqn);
        else
            return false;
    }

    int GetPQN() const { return pqn; }
    OrbitalInfo GetOrbitalInfo() const { return OrbitalInfo(pqn, sym.GetTwoJ(), sym.GetParity()); }

    virtual std::string Name() const override
    {   return GetOrbitalInfo().Name();
    }

    virtual std::string Print() const override
    {   return GetOrbitalInfo().Name();
    }

    virtual void Read(FILE* fp) override
    {
        HamiltonianID::Read(fp);
        pqn = configs->front().begin()->first.PQN();
    }
    
    virtual pHamiltonianID Clone() const
    {   return std::make_shared<SingleOrbitalID>(*this);
    }

protected:
    int pqn;
};

/** NonRelID defines a Hamiltonian generated from RelativisticConfigurations drawn from a single
    NonRelConfiguration.
 */
class NonRelID : public HamiltonianID
{
public:
    NonRelID(): HamiltonianID() {}
    NonRelID(const NonRelConfiguration& config, int two_j, pRelativisticConfigList rconfigs = nullptr):
        HamiltonianID(two_j, config.GetParity(), rconfigs), nrconfig(config) {}
    NonRelID(const NonRelID& other): HamiltonianID(other), nrconfig(other.nrconfig) {}

    virtual bool operator<(const HamiltonianID& other) const override
    {
        if(sym < other.GetSymmetry())
            return true;
        else if(other.GetSymmetry() < sym)
            return false;
        else
        {   const NonRelID* other_nrid = dynamic_cast<const NonRelID*>(&other);
            if(other_nrid)
                return nrconfig < other_nrid->nrconfig;
            else
                return false;
        }
    }

    virtual bool operator==(const HamiltonianID& other) const override
    {
        const NonRelID* other_nrid = dynamic_cast<const NonRelID*>(&other);
        if(other_nrid)
            return (sym == other_nrid->sym) && (nrconfig == other_nrid->nrconfig);
        else
            return false;
    }

    NonRelConfiguration& GetNonRelConfiguration() { return nrconfig; }
    const NonRelConfiguration& GetNonRelConfiguration() const { return nrconfig; }

    virtual std::string Name() const override
    {   return nrconfig.Name() + "." + HamiltonianID::Name();
    }

    virtual std::string Print() const override
    {   return nrconfig.Name() + " " + HamiltonianID::Print();
    }

    virtual void Read(FILE* fp) override
    {   HamiltonianID::Read(fp);
        nrconfig = NonRelConfiguration(configs->front());
    }

    virtual pHamiltonianID Clone() const
    {   return std::make_shared<NonRelID>(*this);
    }

protected:
    NonRelConfiguration nrconfig;
};

#endif
