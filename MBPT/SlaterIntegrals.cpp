#include "Include.h"
#include "SlaterIntegrals.h"

#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"

SlaterIntegrals::SlaterIntegrals(const ExcitedStates* excited_states):
    core(*excited_states->GetCore()), excited(*excited_states), include_valence_sms(false)
{
    if(core.GetNuclearInverseMass())
        include_valence_sms = true;
}

void SlaterIntegrals::UpdateStateIndexes(const ExcitedStates& valence)
{
    if(core.GetNuclearInverseMass())
        include_valence_sms = true;

    NumStates = core.NumStates() + excited.NumStates();
    state_index.clear();
    reverse_state_index.clear();

    ConstStateIterator it_i = core.GetConstStateIterator();
    unsigned int i;

    // Iterate through states, assign in order
    it_i.First(); i = 0;
    while(!it_i.AtEnd())
    {
        core_states.insert(i);
        state_index.insert(std::pair<OrbitalInfo, unsigned int>(OrbitalInfo(it_i.GetState()), i));
        reverse_state_index.insert(std::pair<unsigned int, OrbitalInfo>(i, OrbitalInfo(it_i.GetState())));

        it_i.Next(); i++;
    }

    it_i = excited.GetConstStateIterator();
    it_i.First();
    while(!it_i.AtEnd())
    {
        OrbitalInfo info(it_i.GetState());

        // Check not already present from open shell core
        if(state_index.find(info) == state_index.end())
        {   state_index.insert(std::pair<OrbitalInfo, unsigned int>(info, i));
            reverse_state_index.insert(std::pair<unsigned int, OrbitalInfo>(i, info));
        }

        excited_states.insert(i);
        if(valence.GetState(info))
            valence_states.insert(i);

        it_i.Next(); i++;
    }
}

void SlaterIntegrals::Clear()
{
    state_index.clear();
    reverse_state_index.clear();
    OneElectronIntegrals.clear();
    SMSIntegrals.clear();
    TwoElectronIntegrals.clear();
}

void SlaterIntegrals::Update(const ExcitedStates& valence)
{
    Clear();
    UpdateStateIndexes(valence);

    UpdateOneElectronIntegrals();
    UpdateTwoElectronIntegrals();
}

double SlaterIntegrals::GetOneElectronIntegral(const OrbitalInfo& s1, const OrbitalInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return OneElectronIntegrals.find(i1 * NumStates + i2)->second;
    else
        return OneElectronIntegrals.find(i2 * NumStates + i1)->second;
}

double SlaterIntegrals::GetSMSIntegral(const OrbitalInfo& s1, const OrbitalInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return SMSIntegrals.find(i1 * NumStates + i2)->second;
    else
        return -SMSIntegrals.find(i2 * NumStates + i1)->second;
}

bool SlaterIntegrals::TwoElectronIntegralOrdering(unsigned int& i1, unsigned int& i2, unsigned int& i3, unsigned int& i4) const
{
    bool sms_sign = true;

    // Ordering of indices:
    // (i1 <= i3) && (i2 <= i4) && (i1 <= i2) && (if i1 == i2, then (i3 <= i4))
    // therefore (i1 <= i2 <= i4) and (i1 <= i3)
    if(i3 < i1)
    {   swap(i3, i1);
        sms_sign = !sms_sign;
    }
    if(i4 < i2)
    {   swap(i4, i2);
        sms_sign = !sms_sign;
    }
    if(i2 < i1)
    {   swap(i2, i1);
        swap(i3, i4);
    }
    if((i1 == i2) && (i4 < i3))
        swap(i3, i4);

    return sms_sign;
}

double SlaterIntegrals::GetTwoElectronIntegral(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;
    unsigned int i3 = state_index.find(s3)->second;
    unsigned int i4 = state_index.find(s4)->second;

    bool sms_sign = TwoElectronIntegralOrdering(i1, i2, i3, i4);

    LongKey key = k  * NumStates*NumStates*NumStates*NumStates +
                  i1 * NumStates*NumStates*NumStates +
                  i2 * NumStates*NumStates +
                  i3 * NumStates +
                  i4;

    double radial = 0.;

    if(TwoElectronIntegrals.find(key) != TwoElectronIntegrals.end())
    {
        radial = TwoElectronIntegrals.find(key)->second;
    }
    else
    {   *errstream << "SlaterIntegrals::GetTwoElectronIntegral() failed to find integral."
                   << "\n  key = " << key << "  num_states = " << NumStates 
                   << "\n  R^" << k << " ( " << s1.Name() << " " << s2.Name()
                   << ", " << s3.Name() << " " << s4.Name() << ") :" << std::endl;
    }

    if(include_valence_sms && (k == 1))
    {   double SMS = GetNuclearInverseMass();
        if(SMS)
        {   SMS = SMS * SMSIntegrals.find(i1*NumStates + i3)->second * SMSIntegrals.find(i2*NumStates + i4)->second;
            if(!sms_sign)
                SMS = -SMS;
            radial = radial - SMS;
        }
    }

    return radial;
}
