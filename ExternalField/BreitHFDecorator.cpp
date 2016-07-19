#include "BreitHFDecorator.h"

double BreitHFDecorator::GetMatrixElement(const Orbital& b, const Orbital& a) const
{
    double value = wrapped->GetMatrixElement(b, a);
    if(a.Kappa() != b.Kappa())
        return value;

    bool NON_REL_SCALING = true;

    // Find out whether a is in the core
    const Orbital* current_in_core = dynamic_cast<const Orbital*>(&a);
    if(core->GetState(OrbitalInfo(current_in_core)) == NULL)
        current_in_core = NULL;

    pSpinorFunctionConst p_a(a.shared_from_this());

    // Sum over all core states
    auto cs = core->begin();
    while(cs != core->end())
    {
        pOrbitalConst core_orbital = cs->second;
        double other_occupancy = core->GetOccupancy(OrbitalInfo(core_orbital));

        // Sum over all k
        int k = breit_operator->SetOrbitals(core_orbital, p_a);
        while(k != -1)
        {
            double coefficient = -1./(2. * abs(a.Kappa()));
            coefficient *= MathConstant::Instance()->minus_one_to_the_power((a.TwoJ() - core_orbital->TwoJ())/2 + k);

            // Open shells need to be scaled
            if(other_occupancy != double(2 * abs(core_orbital->Kappa())))
            {
                double ex = 1.;
                if(NON_REL_SCALING)
                {   // Average over non-relativistic configurations
                    if(core_orbital->Kappa() == -1)
                    {
                        if(!current_in_core || (OrbitalInfo(current_in_core) != OrbitalInfo(core_orbital)))
                            ex = other_occupancy/double(2 * abs(core_orbital->Kappa()));
                        else if(k)
                            ex = (other_occupancy - 1.)/double(2 * abs(core_orbital->Kappa()) - 1);
                    }
                    else
                    {   OrbitalInfo pair_info(core_orbital->PQN(), - core_orbital->Kappa() - 1);
                        pOrbitalConst pair_orbital = core->GetState(pair_info);
                        double pair_occupancy = core->GetOccupancy(pair_info);

                        if((!current_in_core && a.L() != core_orbital->L())
                           || (current_in_core && (OrbitalInfo(current_in_core) != OrbitalInfo(core_orbital)) && (OrbitalInfo(current_in_core) != pair_info)))
                            ex = (other_occupancy + pair_occupancy)/double(2 * (abs(core_orbital->Kappa()) + abs(pair_info.Kappa())));
                        else if(k)
                            ex = (other_occupancy + pair_occupancy - 1.)/double(2 * (abs(core_orbital->Kappa()) + abs(pair_info.Kappa())) - 1);
                    }
                }
                else
                {   // Average over relativistic configurations
                    if(!current_in_core || (OrbitalInfo(current_in_core) != OrbitalInfo(core_orbital)))
                        ex = other_occupancy/double(2 * (abs(core_orbital->Kappa())));
                    else if(k)
                        ex = (other_occupancy - 1.)/double(2 * (abs(core_orbital->Kappa())) - 1);
                }

                coefficient = coefficient * ex;
            }

            value += breit_operator->GetMatrixElement(b, *core_orbital) * coefficient;

            k = breit_operator->NextK();
        }

        cs++;
    }
    
    return value;
}

SpinorFunction BreitHFDecorator::CalculateExtraExchange(const SpinorFunction& s) const
{
    bool NON_REL_SCALING = true;

    SpinorFunction exchange(s.Kappa());
    exchange.resize(s.size());

    // Find out whether s is in the core
    const Orbital* current_in_core = dynamic_cast<const Orbital*>(&s);
    if(core->GetState(OrbitalInfo(current_in_core)) == NULL)
        current_in_core = NULL;

    pSpinorFunctionConst p_s(s.shared_from_this());

    // Sum over all core states
    auto cs = core->begin();
    while(cs != core->end())
    {
        pOrbitalConst core_orbital = cs->second;
        double other_occupancy = core->GetOccupancy(OrbitalInfo(core_orbital));

        // Sum over all k
        int k = breit_operator->SetOrbitals(core_orbital, p_s);
        while(k != -1)
        {
            double coefficient = 1./(2. * abs(s.Kappa()));
            coefficient *= MathConstant::Instance()->minus_one_to_the_power((s.TwoJ() - core_orbital->TwoJ())/2 + k);

            // Open shells need to be scaled
            if(other_occupancy != double(2 * abs(core_orbital->Kappa())))
            {
                double ex = 1.;
                if(NON_REL_SCALING)
                {   // Average over non-relativistic configurations
                    if(core_orbital->Kappa() == -1)
                    {
                        if(!current_in_core || (OrbitalInfo(current_in_core) != OrbitalInfo(core_orbital)))
                            ex = other_occupancy/double(2 * abs(core_orbital->Kappa()));
                        else if(k)
                            ex = (other_occupancy - 1.)/double(2 * abs(core_orbital->Kappa()) - 1);
                    }
                    else
                    {   OrbitalInfo pair_info(core_orbital->PQN(), - core_orbital->Kappa() - 1);
                        pOrbitalConst pair_orbital = core->GetState(pair_info);
                        double pair_occupancy = core->GetOccupancy(pair_info);

                        if((!current_in_core && s.L() != core_orbital->L())
                           || (current_in_core && (OrbitalInfo(current_in_core) != OrbitalInfo(core_orbital)) && (OrbitalInfo(current_in_core) != pair_info)))
                            ex = (other_occupancy + pair_occupancy)/double(2 * (abs(core_orbital->Kappa()) + abs(pair_info.Kappa())));
                        else if(k)
                            ex = (other_occupancy + pair_occupancy - 1.)/double(2 * (abs(core_orbital->Kappa()) + abs(pair_info.Kappa())) - 1);
                    }
                }
                else
                {   // Average over relativistic configurations
                    if(!current_in_core || (OrbitalInfo(current_in_core) != OrbitalInfo(core_orbital)))
                        ex = other_occupancy/double(2 * (abs(core_orbital->Kappa())));
                    else if(k)
                        ex = (other_occupancy - 1.)/double(2 * (abs(core_orbital->Kappa())) - 1);
                }

                coefficient = coefficient * ex;
            }

            exchange += breit_operator->ApplyTo(*core_orbital, s.Kappa()) * coefficient;

            k = breit_operator->NextK();
        }

        cs++;
    }
    
    return exchange;
}
