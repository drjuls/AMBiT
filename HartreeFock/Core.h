#ifndef CORE_H
#define CORE_H

#include "StateManager.h"
#include "DiscreteStateIterator.h"
#include "Atom/Debug.h"
#include "MBPT/SigmaPotential.h"
#include "ContinuumState.h"

class Core : public StateManager
{
    /** The core is closed; all electron shells are filled.
        Core does all of the Hartree-Fock calculation.
     */
public:
    friend class DiscreteStateIterator;
    friend class ConstDiscreteStateIterator;

public:     // Methods for StateManager role
    Core(Lattice* lat, unsigned int atomic_number, int ion_charge);
    virtual ~Core() {}

    /** This method must be called before any other, except for Read() or GetDebugOptions(). */
    virtual void Initialise();

    /** These override the methods in StateManager to return a DiscreteState. */
    virtual DiscreteState* GetState(const StateInfo& info);
    virtual const DiscreteState* GetState(const StateInfo& info) const;
    
    virtual DiscreteStateIterator GetDiscreteStateIterator();
    virtual ConstDiscreteStateIterator GetConstDiscreteStateIterator() const;

    virtual void Write(FILE* fp) const;

    /** Use instead of initialisation to save computational time. */
    virtual void Read(FILE* fp);

public:     // Methods for setting physical parameters of core
    /** Get the number of protons in the nucleus (atomic number) */
    double GetZ() const { return Z; }

    /** Get the charge of the core (ionisation degree) */
    double GetCharge() const { return Charge; }

    /** Nuclear parameters */
    void SetNuclearRadius(double rnuc);
    void SetNuclearThickness(double dnuc);
    void SetNuclearInverseMass(double inv_mass);
    void SetVolumeShiftParameter(double vol_shift);

    double GetNuclearRadius() const { return NuclearRadius; }
    double GetNuclearThickness() const { return NuclearThickness; }
    double GetNuclearInverseMass() const { return NuclearInverseMass; }
    double GetVolumeShiftParameter() const { return VolumeShiftParameter; }

    /** Create additional volume shift potential. This is added to the nuclear potential
        with a scaling factor (VolumeShiftParameter).
        Supply the difference in nuclear radii, dR, to obtain for nuclear potential
            dV(R, r)   V(R + dR, r) - V(R, r)
            -------- = ----------------------
               dR                dR
      */
    void CalculateVolumeShiftPotential(double radius_difference);

public:     // Methods for Hartree-Fock calculations and potentials
    /** Set energy tolerances */
    void SetEnergyTolerance(double tolerance);
    double GetEnergyTolerance() { return StateParameters::EnergyTolerance; }

    /** Do self consistent Hartree Fock proceedure until convergency is reached. */
    void Update();
    /** If convergency is not reached, return false. */
    bool UpdateGreens();

    /** The HF Potential is the direct potential of the core including open shells. */
    virtual std::vector<double> GetHFPotential() const;
    virtual std::vector<double> GetLocalExchangeApproximation() const;

    /** Calculate the exchange part of the interaction between the current state and all
        of the core states. "exchange" is resized to current.size().
     */
    void CalculateExchange(const State& current, CoupledFunction& exchange, const SigmaPotential* sigma = NULL, double sigma_amount = 1.) const;

    /** Calculate a new excited state in the closed core potential. */
    virtual unsigned int CalculateExcitedState(State* s) const;

    /** Update an existing excited state, in case of changed core or addition of sigma. */
    unsigned int UpdateExcitedState(State* s, const SigmaPotential* sigma = NULL, double sigma_amount = 1.) const;

public:
    /** Methods for open shell core. The core calculates states in the V^n scheme.
        Can also toggle modes for V^(n-1) and V^(n-c), which is a closed shell core.
      */

    /** Remove one electron from outer shell. Good for making a V^(n-1) core. */
    virtual void Ionise(StateInfo removed_electron);
    
    /** Remove all electrons until just the closed shell core is left. */
    virtual void ToggleClosedShellCore();

    /** Put all the old electrons back onto the core to get the V^n core. */
    virtual void ToggleOpenShellCore();

    /** Test whether a state is in the open shells of the core. */
    virtual bool IsOpenShellState(const StateInfo& info) const;

protected:
    /** Iterate an existing state in the current HF potential until the energy converges.
        It includes the exchange interaction with all core electrons,
            scaled by exchange_amount (which should be between 0 and 1, inclusive).
        If exchange_amount < 0, then use the local approximation to exchange for the
            first iteration, then the full exchange thereafter.
        Returns the number of iterations necessary to achieve convergence.
        If convergency is not reached, it returns a value >= StateParameters::MaxHFIterations.
     */
    unsigned int CalculateDiscreteState(DiscreteState* s, double exchange_amount = 1., const SigmaPotential* sigma = NULL, double sigma_amount = 1.) const;

    /** Solve the Dirac eqn in the current HF potential once and adjust the energy
        to renormalise. The exchange may be resized to the new size of s.
        Returns the change in the energy needed for renormalisation (delta_E).
     */
    double IterateDiscreteStateGreens(DiscreteState* s, CoupledFunction* exchange) const;

    /** Calculate a new continuum state in the HF Potential. */
    unsigned int CalculateContinuumState(ContinuumState* s) const;

    /** Update the Hartree-Fock potential, which includes the potential due to the nucleus,
        as well as the direct part of all of the core electrons.
        The first_build parameter is a hack, due to the need to treat the LocalExchangeApproximation
        slightly differently on the first build (before starting HF iterations). When this is true,
        the LocalExchangeApproximation is added to HFPotential in this function.
     */
    void UpdateHFPotential(double proportion_new = 1., bool first_build = false);

    void UpdateHFPotentialNoSelfInteraction(const StateInfo& current, double proportion_new = 1.);
    void CalculateExchangeNoSelfInteraction(const State& current, CoupledFunction& exchange) const;

    /** Simply extend the HF potential (and the closed core potential) to the size of the lattice. */
    void ExtendPotential() const;

    /** Update the coulomb potential due to the nucleus, inside the nucleus */
    void UpdateNuclearPotential();

    /** Calculate the nuclear density. */
    std::vector<double> CalculateNuclearDensity(double radius, double thickness) const;

    /** Build the first approximation to the atom state */
    void BuildFirstApproximation();

    /** Used by BuildFirstApproximation to get the configuration, if it has been done before. */
    bool GetConfigData();

    /** Used by BuildFirstApproximation to get states which don't fit into the standard
        electron filling scheme. Returns number of electrons allocated.
     */
    unsigned int GetSpecialStates();

protected:
    double NuclearRadius;
    double NuclearThickness;
    double NuclearInverseMass;
    std::vector<double> NuclearPotential;

    /** The volume shift parameter arbitrarily increases the size of the volume shift. */
    double VolumeShiftParameter;
    std::vector<double> VolumeShiftPotential;

    /** HFPotential is mutable in case it's size needs to be increased. */
    mutable std::vector<double> HFPotential;
    mutable std::vector<double> LocalExchangeApproximation;

    StateSet OpenShellStorage;  // Store currently unused open shell states
    std::map<StateInfo, double> OpenShellStates; // Original occupancies

public:
    class StateParameters
    {
    public:
        /** Maximum number of Hartree Fock iterations before admitting that
            there is a convergency problem
         */
        const static unsigned int MaxHFIterations;
        /** Cutoff wavefunction when the ratio
              f[last_point]/f_max < WavefunctionTolerance
         */
        static double WavefunctionTolerance;
        /** Accuracy of Hartree Fock energy
         */
        static double EnergyTolerance;
    };
};

inline void Core::SetNuclearRadius(double rnuc)
{
    NuclearRadius = rnuc;
    UpdateNuclearPotential();
    Update();
}

inline void Core::SetNuclearThickness(double dnuc)
{   
    NuclearThickness = dnuc;
    UpdateNuclearPotential();
    Update();
}

inline void Core::SetNuclearInverseMass(double inv_mass)
{
    NuclearInverseMass = inv_mass;
}

inline void Core::SetVolumeShiftParameter(double vol_shift)
{
    VolumeShiftParameter = vol_shift;
    //if(VolumeShiftPotential.size())
    //    Update();
}

inline void Core::SetEnergyTolerance(double tolerance)
{
    StateParameters::EnergyTolerance = tolerance;
    Update();
}

#endif
