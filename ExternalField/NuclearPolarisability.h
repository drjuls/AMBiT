#ifndef NUCLEAR_POLARISABILITY_H
#define NUCLEAR_POLARISABILITY_H

#include "HartreeFock/HFOperator.h"
#include "HartreeFock/LocalPotentialDecorator.h"

namespace Ambit
{
/** The nuclear polarisability has the form
    \f[
        V(r) = \frac{-\alpha_E}{2 (r^4 + (\frac{\alpha \pi^2}{\sqrt{2}R})^4)}
    \f]
    and
    \f[
        R = \frac{19}{6} + 5 \log(\frac{2\bar{E}}{m_e c^2})
    \f]
    where \f$\alpha_E\f$ is the nuclear electric polarisability and \f$\bar{E}\f$ is the mean nuclear excitation energy.
 */
class NuclearPolarisability: public HFOperatorDecorator<LocalPotentialDecorator, NuclearPolarisability>
{
public:
    NuclearPolarisability(pHFOperator wrapped_hf, double alphaE, double nuclear_excitation_energy_MeV, pIntegrator integration_strategy = pIntegrator());

protected:
    void GeneratePotential();
    double Ebar;
};

}
#endif
