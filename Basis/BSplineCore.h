#ifndef BSPLINE_CORE_H
#define BSPLINE_CORE_H

#include "BSplineGrid.h"
#include "HartreeFock/Core.h"

class BSplineCore
{
public:
    BSplineCore(BSplineGrid* lat, const Core* core);
    virtual ~BSplineCore(void) {}

    /** Calculate exchange of BSpline (numerated according to grid parameters) with the core.
     *  The spline could be either an upper or lower component (f or g)
     *  which affects the density function.
     */
    void CalculateExchange(unsigned int bspline, int kappa, bool upper, CoupledFunction& exchange);

    /** Get direct HF potential on the bspline grid. */
    std::vector<double> GetPotential() { return potential; }

protected:
    void Read(FILE* fp) {}
    void Write(FILE* fp) const {}

protected:
    BSplineGrid* grid;
    const Core* hfcore;
    std::vector<double> potential;

    unsigned int order;     // interpolation order
};

#endif
