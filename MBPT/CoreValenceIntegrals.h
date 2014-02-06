#ifndef CORE_VALENCE_INTEGRALS_H
#define CORE_VALENCE_INTEGRALS_H

#include "SlaterIntegrals.h"

class CoreValenceIntegrals : public SlaterIntegrals
{
    /** Class to hold Slater integrals for use in MBPT calculation. */
public:
    CoreValenceIntegrals(pExcitedStatesConst excited_states):
        SlaterIntegrals(excited_states)
    {}
    virtual ~CoreValenceIntegrals() {}

    /** Calculate number of one-electron and two-electron integrals that will be stored.
        Return total.
     */
    virtual unsigned int GetStorageSize(const ExcitedStates& valence);

protected:
    virtual void UpdateOneElectronIntegrals();
    virtual void UpdateTwoElectronIntegrals();
};

#endif
