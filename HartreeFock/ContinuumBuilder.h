#ifndef CONTINUUM_BUILDER_H
#define CONTINUUM_BUILDER_H

#include "Core.h"

#ifndef CONTINUUM_NORMALISATION_ENUM
#define CONTINUUM_NORMALISATION_ENUM
    enum ContinuumNormalisation { Flambaum, Cowan, Unitary };
    /* Amplitudes as r -> Infinity
        Flambaum normalization:  A = 2 * Pi^(-1/2) * E^(1/2)
        Cowan normalization:     A = Pi^(-1/2) * (2/E)^(1/4)
        Unitary normalization:   A = 1
     */
#endif

class ContinuumBuilder
{
    /** Class to create continuum wavefunctions.
        ContinuumBuilder lets us modify the Hartree-Fock core, or even create
        an entirely new one, if required.
        The lattice may be extended or modified at will to help create continuum
        wavefunctions, and the resulting wavefunctions are interpolated onto the
        standard lattice.
     */
public:
    ContinuumBuilder(): lattice(NULL), core(NULL), norm_type(Cowan) {}
    /** Copy lattice and core from other_core. */
    ContinuumBuilder(const Core* other_core);
    virtual ~ContinuumBuilder();

    void CopyLattice(const Lattice* lat);

    /** Create a new regular Lattice. These work best because they tend towards
        even spacing as r->Infinity.
        For other possibilities create externally and use CopyLattice. 
     */
    void CreateNewLattice(unsigned int numpoints, double r_min, double r_max);

    /** Copy existing HF core.
        If import_lattice is true, function replaces lattice,
            otherwise interpolates core onto current lattice.
     */
    void CopyCore(const Core* other_core, bool import_lattice = false);

    /** Use GetCore() to make modifications (such as ionization of an electron). */
    Core* GetCore() { return core; }
    const Core* GetCore() const { return core; }

public:
    /** Set normalisation of continuum wavefunctions. */
    void SetNormalisationType(ContinuumNormalisation norm)
    {   norm_type = norm;
    }
    
    ContinuumNormalisation GetNormalisationType() const
    {   return norm_type;
    }

    /** Create a new continuum wavefunction in the field of the core.
        Optionally, interpolate continuum wavefunction onto external_lattice.
        PRE: core and lattice must already exist.
        Return number of loops required for HF convergence.
        Return value of 0 means start_sine was not reached and method failed.
     */
    unsigned int CalculateContinuumWave(pContinuumWave s, Lattice* external_lattice = NULL) const;

    /** Read a continuum wavefunction from file and interpolate onto external_lattice.
        The normalisation type of the continuum builder is set to Unitary.
        Return success.
     */
    bool ReadContinuumWave(pContinuumWave s, Lattice* external_lattice, const std::string& upper_file, const std::string& lower_file);

protected:
    Lattice* lattice;
    Core* core;

    ContinuumNormalisation norm_type;
};

#endif
