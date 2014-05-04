#ifndef EXCITED_STATES_H
#define EXCITED_STATES_H

#include "HartreeFock/OrbitalMap.h"
#include "HartreeFock/Core.h"
#include "HartreeFock/SigmaPotential.h"
#include <set>
#include <map>
#include <vector>

typedef std::map<int, SigmaPotential*> SigmaMap;
typedef std::map<OrbitalInfo, double> SigmaAmount;

class CoreMBPTCalculator;

class ExcitedStates : public OrbitalMap
{
public:
    ExcitedStates(pLattice lattice);
    virtual ~ExcitedStates();

public:     // Methods for managing sigma
    /** Set an identifier to use in the sigma file names. */
    void SetIdentifier(const std::string& s) { identifier = s; }
    
    /** Clear all sigmas. */
    void ClearSigmas();

    /** Calculate the matrix element  < s | Sigma | s >. */
    double GetSigmaMatrixElement(const OrbitalInfo& info) const;

    /** Get a copy of the state, including iterating sigma (if available). */
    virtual Orbital GetStateWithSigma(const OrbitalInfo& info) const;

    /** Create a sigma operator for the given state to second order.
        This operator will henceforth be included in the exchange potential of
        the state (to turn this off set SigmaAmount to zero).
        Return value is the energy of the state including second order matrix elements.
     */
    double CreateSecondOrderSigma(const OrbitalInfo& info, const CoreMBPTCalculator& mbpt);

    /** Retreive a sigma operator that already exists. Return true if successful. */
    bool RetrieveSecondOrderSigma(const OrbitalInfo& info);

    /** Set SigmaAmount to give correct energy (in Hartree units).
     */
    void SetEnergyViaSigma(const OrbitalInfo& info, double energy);

    /** Get/Set the amount of sigma to be used when calculating state for which sigma exists. */
    void SetSigmaAmount(const OrbitalInfo& info, double amount);
    double GetSigmaAmount(const OrbitalInfo& info) const;

    /** Test for orthogonality of orbitals with each other and the core.
        Return largest overlap.
     */
    double TestOrthogonality(pCoreConst core) const;

protected:
    /** current.f = r * previous.f
        PRE: previous.kappa = current.kappa
     */
    void MultiplyByR(pOrbitalConst previous, pOrbital current) const;

    /** current.f = sin(kr) * previous.f, where k = Pi/R_max
        PRE: previous.kappa = current.kappa
     */
    void MultiplyBySinR(pOrbitalConst previous, pOrbital current) const;

    /** current.f = r * sin(kr) * previous.f
        PRE: current.l = previous.l + 1
     */
    void MultiplyByRSinR(pOrbitalConst previous, pOrbital current) const;

    /** Orthogonalise to all states that have the same kappa and principal quantum number
        less than current (both in core and in excited states).
        PRE: all states with smaller pqn must already be orthogonal
     */
    void Orthogonalise(pOrbital current) const;
    
    /** Delete all currently stored states and sigma potentials. */
    virtual void Clear();

protected:
    std::string identifier;  // needed to store sigmas
    SigmaMap SecondOrderSigma;
    SigmaAmount SecondOrderAmount;
};

typedef boost::shared_ptr<ExcitedStates> pExcitedStates;
typedef boost::shared_ptr<const ExcitedStates> pExcitedStatesConst;

#endif
