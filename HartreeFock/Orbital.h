#ifndef ORBITAL_H
#define ORBITAL_H

#include "Universal/Lattice.h"
#include "Universal/SpinorFunction.h"
#include "OpIntegrator.h"
#include <memory>
#include <string>
#include <math.h>
#include <stdlib.h>

class OrbitalInfo;

/** Orbital is for solutions of a single-particle Hamiltonian.
    That is, it is a SpinorFunction with an energy eigenvalue E:
        \f[ h|\psi> = E|\psi> \f]
    It also contains a principal quantum number, which may or may not correspond
    to the spectroscopic pqn.
    Some of the functions assume that the orbital is explicitly bounded and hence normalisable, integrable, etc.
 */
class Orbital : public SpinorFunction
{
public:
    Orbital(const OrbitalInfo& info);
    Orbital(int kappa, int pqn = 0, double energy = 0.0, unsigned int size = 0);
    Orbital(const Orbital& other);
    Orbital(Orbital&& other);
    virtual ~Orbital() {}
    
    const Orbital& operator=(const Orbital& other);
    Orbital& operator=(Orbital&& other);
    
    virtual double Energy() const;
    
    /** Nu is the effective principal quantum number, E = -1/(2 nu^2). */
    virtual double Nu() const;
    virtual int PQN() const;
    
    virtual void SetEnergy(double Energy);
    virtual void SetNu(double Nu);
    virtual void SetPQN(int PQN);
    
    virtual std::string Name() const;
    
    const Orbital& operator*=(double scale_factor);
    Orbital operator*(double scale_factor) const;
    const Orbital& operator+=(const Orbital& other);
    const Orbital& operator-=(const Orbital& other);
    Orbital operator+(const Orbital& other) const;
    Orbital operator-(const Orbital& other) const;
    const Orbital& operator*=(const RadialFunction& chi);
    Orbital operator*(const RadialFunction& chi) const;
    
    /** Store the state. File pointer fp must be open and writable. */
    virtual void Write(FILE* fp) const;
    
    /** Read the state from file. fp must be open and readable. */
    virtual void Read(FILE* fp);
    
    /** Print state to file, optionally printing lattice. Return success. */
    bool Print(pLattice lattice = pLattice()) const;
    bool Print(const std::string& filename, pLattice lattice = pLattice()) const;
    bool Print(FILE* fp, pLattice lattice = pLattice()) const;

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
    double Norm(pOPIntegrator integrator) const;
    
    /** Scale the orbital so that it is normalised to "norm". */
    void ReNormalise(pOPIntegrator integrator, double norm = 1.);
    
    /** Count the number of nodes of the wavefunction. */
    unsigned int NumNodes() const;

protected:
    int pqn;
    double energy;
};

typedef std::shared_ptr<Orbital> pOrbital;
typedef std::shared_ptr<const Orbital> pOrbitalConst;

typedef Orbital ContinuumWave;
typedef pOrbital pContinuumWave;
typedef pOrbitalConst pContinuumWaveConst;

#endif

