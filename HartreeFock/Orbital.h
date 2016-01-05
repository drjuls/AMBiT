#ifndef ORBITAL_H
#define ORBITAL_H

#include "Universal/Lattice.h"
#include "Universal/SpinorFunction.h"
#include "OpIntegrator.h"
#include <memory>
#include <string>
#include <math.h>
#include <stdlib.h>

class Orbital;
class OrbitalInfo;

/** Orbital is for solutions of a single-particle Hamiltonian.
    That is, it is a SpinorFunction with an energy eigenvalue E:
        \f[ h|\psi> = E|\psi> \f]
    It also contains a principal quantum number, which may or may not correspond
    to the spectroscopic pqn.
    Some of the functions assume that the orbital is explicitly bounded and hence normalisable, integrable, etc.
 */
class OrbitalBase : public SpinorFunction
{
public:
    OrbitalBase(const OrbitalInfo& info);
    OrbitalBase(int kappa, int pqn = 0, double energy = 0.0, unsigned int size = 0);
    virtual ~OrbitalBase() {}

    virtual std::shared_ptr<Orbital> Clone() const = 0; // { return std::make_shared<OrbitalBase>(*this); };

    virtual double Energy() const;
    
    /** Nu is the effective principal quantum number, E = -1/(2 nu^2). */
    virtual double Nu() const;
    virtual int PQN() const;
    
    virtual void SetEnergy(double Energy);
    virtual void SetNu(double Nu);
    virtual void SetPQN(int PQN);
    
    virtual std::string Name() const;
    
//    const OrbitalBase& operator*=(double scale_factor);
//    OrbitalBase operator*(double scale_factor) const;
//    const OrbitalBase& operator+=(const OrbitalBase& other);
//    const OrbitalBase& operator-=(const OrbitalBase& other);
//    OrbitalBase operator+(const OrbitalBase& other) const;
//    OrbitalBase operator-(const OrbitalBase& other) const;
//    const OrbitalBase& operator*=(const RadialFunction& chi);
//    OrbitalBase operator*(const RadialFunction& chi) const;

    /** Store the state. File pointer fp must be open and writable. */
    virtual void Write(FILE* fp) const;
    
    /** Read the state from file. fp must be open and readable. */
    virtual void Read(FILE* fp);
    
    /** Print state to file, optionally printing lattice. Return success. */
    virtual bool Print(pLattice lattice = pLattice()) const;
    virtual bool Print(const std::string& filename, pLattice lattice = pLattice()) const;
    virtual bool Print(FILE* fp, pLattice lattice = pLattice()) const;

public: // These functions assume that the orbital is explicitly bounded
    /** Check that the ratio
          f[size()-1]/f_max < tolerance
        and
          f[size()-2]/f_max > tolerance
        and if they are not, enlarge or reduce the size of the wavefunctions.
        Return true if the size was correct, false if there was a change.
     */
    bool CheckSize(pLattice lattice, double tolerance);
    
    /** Get current normalisation of the orbital. */
    virtual double Norm(pOPIntegrator integrator) const;
    
    /** Scale the orbital so that it is normalised to "norm". */
    virtual void ReNormalise(pOPIntegrator integrator, double norm = 1.);
    
    /** Count the number of nodes of the wavefunction. */
    unsigned int NumNodes() const;

protected:
    int pqn;
    double energy;
};

//typedef std::shared_ptr<OrbitalBase> pOrbital;
//typedef std::shared_ptr<const OrbitalBase> pOrbitalConst;

template <typename Base, typename Derived>
class OrbitalTemplate : public Base
{
public:
    using Base::Base;
    OrbitalTemplate(const Base& other): Base(other) {}

    virtual std::shared_ptr<Orbital> Clone() const override;

    /** Multiply all points of all vectors (f, g, df, dg) by the scale factor. */
    const Derived& operator*=(double scale_factor);
    Derived operator*(double scale_factor) const;

    /** Adding or subtracting two spinor functions can only occur if both have same angular part. */
    const Derived& operator+=(const Derived& other);
    const Derived& operator-=(const Derived& other);
    Derived operator+(const Derived& other) const;
    Derived operator-(const Derived& other) const;

    /** Multiply spinor function by another radial function (assumed zero outside range). */
    const Derived& operator*=(const RadialFunction& chi);
    Derived operator*(const RadialFunction& chi) const;

protected:
    typedef OrbitalTemplate<Base, Derived> BaseClass;
};

class Orbital : public OrbitalTemplate<OrbitalBase, Orbital>
{
public:
    using BaseClass::BaseClass;
};

typedef std::shared_ptr<Orbital> pOrbital;
typedef std::shared_ptr<const Orbital> pOrbitalConst;

typedef Orbital ContinuumWave;
typedef std::shared_ptr<Orbital> pContinuumWave;
typedef std::shared_ptr<const Orbital> pContinuumWaveConst;

// Template functions below

template <typename Base, typename Derived>
pOrbital OrbitalTemplate<Base, Derived>::Clone() const
{
    std::shared_ptr<Derived> ret = std::make_shared<Derived>(static_cast<Derived const &>(*this));
    return std::static_pointer_cast<Orbital>(ret);
}

template <typename Base, typename Derived>
auto OrbitalTemplate<Base, Derived>::operator*=(double scale_factor) -> const Derived&
{
    if(scale_factor != 1.)
    {
        for(auto it = this->f.begin(); it != this->f.end(); it++)
        {   (*it) *= scale_factor;
        }
        for(auto it = this->g.begin(); it != this->g.end(); it++)
        {   (*it) *= scale_factor;
        }
        for(auto it = this->dfdr.begin(); it != this->dfdr.end(); it++)
        {   (*it) *= scale_factor;
        }
        for(auto it = this->dgdr.begin(); it != this->dgdr.end(); it++)
        {   (*it) *= scale_factor;
        }
    }

    return static_cast<const Derived&>(*this);
}

template <typename Base, typename Derived>
auto OrbitalTemplate<Base, Derived>::operator*(double scale_factor) const -> Derived
{
    Derived ret(static_cast<const Derived&>(*this));
    return (ret *= scale_factor);
}

template <typename Base, typename Derived>
auto OrbitalTemplate<Base, Derived>::operator+=(const Derived& other) -> const Derived&
{
    if(this->size() < other.size())
        this->resize(other.size());

        for(unsigned int i = 0; i < other.size(); i++)
        {   this->f[i] += other.f[i];
            this->g[i] += other.g[i];
            this->dfdr[i] += other.dfdr[i];
            this->dgdr[i] += other.dgdr[i];
        }

    return static_cast<const Derived&>(*this);
}

template <typename Base, typename Derived>
auto OrbitalTemplate<Base, Derived>::operator-=(const Derived& other) -> const Derived&
{
    (*this) += other * (-1.0);
    return static_cast<const Derived&>(*this);
}

template <typename Base, typename Derived>
auto OrbitalTemplate<Base, Derived>::operator+(const Derived& other) const -> Derived
{
    Derived ret(*this);
    return ret += other;
}

template <typename Base, typename Derived>
auto OrbitalTemplate<Base, Derived>::operator-(const Derived& other) const -> Derived
{
    Derived ret(*this);
    return ret -= other;
}

template <typename Base, typename Derived>
auto OrbitalTemplate<Base, Derived>::operator*=(const RadialFunction& chi) -> const Derived&
{
    // Outside range of chi, chi is assumed to be zero.
    if(chi.size() < this->size())
        this->resize(chi.size());

        for(unsigned int i = 0; i < this->size(); i++)
        {   this->f[i] *= chi.f[i];
            this->g[i] *= chi.f[i];
            this->dfdr[i] = this->f[i] * chi.dfdr[i] + this->dfdr[i] * chi.f[i];
            this->dgdr[i] = this->g[i] * chi.dfdr[i] + this->dgdr[i] * chi.f[i];
        }

    return static_cast<const Derived&>(*this);
}

template <typename Base, typename Derived>
auto OrbitalTemplate<Base, Derived>::operator*(const RadialFunction& chi) const -> Derived
{
    Derived ret(*this);
    return (ret *= chi);
}

#endif

