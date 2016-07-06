#include "SpecificMassShiftDecorator.h"
#include "TwoBodySMSOperator.h"
#include "Include.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"

SpecificMassShiftDecorator::SpecificMassShiftDecorator(pHFOperator wrapped_hf, bool nonrel, bool nonrel_include_lower):
    BaseDecorator(wrapped_hf), sms_operator(nullptr)
{
    pHartreeY zero = std::make_shared<HartreeYBase>();
    if(nonrel)
    {
        sms_operator = std::make_shared<TwoBodySMSOperator>(zero, nonrel_include_lower, integrator);
    }
    else
    {
        double Zalpha = Z * physicalConstant->GetAlpha();
        sms_operator = std::make_shared<TwoBodySMSOperator>(zero, Zalpha, integrator);
    }

    sms_operator->SetInverseMass(0.0);
}

SpinorFunction SpecificMassShiftDecorator::CalculateExtraExchange(const SpinorFunction& s) const
{
    bool NON_REL_SCALING = true;

    SpinorFunction exchange(s.Kappa());

    if(sms_operator->GetInverseMass() == 0)
        return exchange;

    // Find out whether s is in the core
    const Orbital* current_in_core = dynamic_cast<const Orbital*>(&s);
    if(core->GetState(OrbitalInfo(current_in_core)) == nullptr)
        current_in_core = NULL;

    pSpinorFunctionConst p_s(s.shared_from_this());

    // Sum over all core states
    auto cs = core->begin();
    while(cs != core->end())
    {
        pOrbitalConst core_orbital = cs->second;
        double other_occupancy = core->GetOccupancy(OrbitalInfo(core_orbital));

        int k = sms_operator->SetOrbitals(p_s, core_orbital);

        // Sum over all k
        while(k != -1)
        {
            double coefficient = MathConstant::Instance()->Electron3j(s.TwoJ(), core_orbital->TwoJ(), k);
            coefficient = (2 * abs(core_orbital->Kappa())) * coefficient * coefficient;

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

            exchange -= sms_operator->ApplyTo(*core_orbital, s.Kappa()) * coefficient;

            k = sms_operator->NextK();
        }
        cs++;
    }

    return exchange;
}
