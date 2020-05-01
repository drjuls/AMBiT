#ifndef CUSTOM_BASIS_H
#define CUSTOM_BASIS_H

#include "HartreeFock/HFOperator.h"
#include "HartreeFock/NonRelInfo.h"

namespace Ambit
{
class CustomBasis
{
public:
    CustomBasis(pHFOperator hf_core);
    virtual ~CustomBasis() {}

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

protected:
    pHFOperator hf;
    RadialFunction potential;
};

}
#endif
