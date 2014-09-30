#ifndef NUCLEUS_DECORATOR_H
#define NUCLEUS_DECORATOR_H

#include "LocalPotentialDecorator.h"

/** Add nuclear field to the direct potential. Assumed Fermi disribution of proton density:
     \f[ \rho(r) = \frac{1}{Exp[B(r - R)] + 1}, \
               B = 4ln(3)/t
     \f]
    where R is the nuclear radius and t a thickness parameter (default 2.3).
    From this a potential is created that is then normalised to Z.
 */
class NucleusDecorator: public LocalPotentialDecorator
{
public:
    NucleusDecorator(pHFOperator wrapped_hf, pOPIntegrator integration_strategy = pOPIntegrator());

    /** Set parameters of proton distribution (in Fermi). */
    virtual void SetFermiParameters(double radius, double thickness = 2.3);
    virtual double GetNuclearRadius() const    { return nuclear_radius;    }    //!< Radius of nucleus in Fermi.
    virtual double GetNuclearThickness() const { return nuclear_thickness; }    //!< Nuclear thickness in Fermi.

protected:
    virtual RadialFunction CalculateNuclearDensity(double radius, double thickness) const;

protected:
    double nuclear_radius;
    double nuclear_thickness;
};

typedef std::shared_ptr<NucleusDecorator> pNucleusDecorator;
typedef std::shared_ptr<const NucleusDecorator> pNucleusDecoratorConst;

#endif
