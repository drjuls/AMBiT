#ifndef CORE_H
#define CORE_H

#include "StateManager.h"
#include "StateIterator.h"
#include "Atom/Debug.h"
#include "SigmaPotential.h"
#include "ContinuumState.h"

class Core : public StateManager
{
    /** The core is closed; all electron shells are filled.
        Core does all of the Hartree-Fock calculation.
     */
public:     // Methods for StateManager role
    Core(Lattice* lat, unsigned int atomic_number, int ion_charge);
    /** Copy all states, interpolating onto new_lattice (if supplied and required).
        If new_lattice is NULL, then lattice = other.lattice (same object, no copy made).
     */
    Core(const Core& other, Lattice* new_lattice = NULL);
    virtual ~Core() {}

    /** This method must be called before any other, except for Read() or GetDebugOptions(). */
    virtual void Initialise(std::string configuration = "");

    /** Writes stored open shell states also. */
    virtual void Write(FILE* fp) const;

    /** This redefines StateManager::Read() because it clears the existing core states
        before reading in all of the stored states. The number of core states is constant
        for any given ion configuration.
        Can use this instead of initialisation to save computational time.
     */
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

    /** We can approximate valence-core correlation effects using polarisability:
           U = (-1/2) * Polarisability/(r^4 + a^4)
        where a is the radius of the core. This models the effect of the electric
        field of the valence electron on the core.
     */
    void SetPolarisability(double pol);
    double GetPolarisability() const { return Polarisability; }

    void CalculateClosedShellRadius();
    double GetClosedShellRadius() const { return ClosedShellRadius; }

public:     // Methods for Hartree-Fock calculations and potentials
    /** Set energy tolerances */
    void SetEnergyTolerance(double tolerance);
    double GetEnergyTolerance() const { return StateParameters::EnergyTolerance; }

    /** Do self consistent Hartree Fock proceedure until convergency is reached. */
    void Update();

    /** Simply extend the HF potential (and the closed core potential) to the size of the lattice. */
    void ExtendPotential() const;

    /** The HF Potential is the direct potential of the core including open shells. */
    virtual std::vector<double> GetHFPotential() const;
    virtual const std::vector<double>& GetConstHFPotential() const;
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

    /** Test whether the core has open shells. */
    virtual bool IsOpenShellCore() const;

    /** Set open shell state occupancy. */
    virtual void SetOpenShellState(const StateInfo& info, double occupancy);

protected:
    /** Delete all currently stored states. */
    virtual void Clear();

    /** Iterate an existing state in an approximate potential until the energy converges.
        The potential is direct + local exchange approximation.
        Returns the number of iterations necessary to achieve convergence.
        If convergency is not reached, it returns a value >= StateParameters::MaxHFIterations.
     */
    unsigned int ConvergeStateApproximation(DiscreteState* s, bool include_exch = true) const;

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

    /** Update the coulomb potential due to the nucleus, inside the nucleus */
    void UpdateNuclearPotential();

    /** Calculate the nuclear density. */
    std::vector<double> CalculateNuclearDensity(double radius, double thickness) const;

    /** Build the first approximation to the atom state */
    void BuildFirstApproximation(std::string configuration);

    /** Used by BuildFirstApproximation to get the configuration, if it has been done before. */
    std::string GetConfigData();

protected:
    double NuclearRadius;
    double NuclearThickness;
    double NuclearInverseMass;
    std::vector<double> NuclearPotential;

    /** The volume shift parameter arbitrarily increases the size of the volume shift. */
    double VolumeShiftParameter;
    std::vector<double> VolumeShiftPotential;

    /** Core polarisability, defined by dE = alpha * E^2. */
    double Polarisability;

    /** Radius of the closed shell core */
    double ClosedShellRadius;

    /** HFPotential is mutable in case it's size needs to be increased. */
    mutable std::vector<double> HFPotential;
    mutable std::vector<double> LocalExchangeApproximation;

    StateSet OpenShellStorage;  // Store currently unused open shell states.
                                // Only storing each state in one place helps memory management.
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
        static double FirstBuildEnergyTolerance;
    };
};

inline void Core::SetNuclearRadius(double rnuc)
{
    NuclearRadius = rnuc;
}

inline void Core::SetNuclearThickness(double dnuc)
{   
    NuclearThickness = dnuc;
}

inline void Core::SetNuclearInverseMass(double inv_mass)
{
    NuclearInverseMass = inv_mass;
}

inline void Core::SetVolumeShiftParameter(double vol_shift)
{
    VolumeShiftParameter = vol_shift;
}

inline void Core::SetPolarisability(double pol)
{
    Polarisability = pol;
}

inline void Core::SetEnergyTolerance(double tolerance)
{
    StateParameters::EnergyTolerance = tolerance;
}

#endif
