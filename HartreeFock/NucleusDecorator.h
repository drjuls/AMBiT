#ifndef NUCLEUS_DECORATOR_H
#define NUCLEUS_DECORATOR_H

#include "LocalPotentialDecorator.h"

/** Add nuclear field to the direct potential. Assumed Fermi disribution of proton density:
     \f[ \rho(r) = \frac{1}{Exp[B(r - R)] + 1}, \
               B = 4ln(3)/t
     \f]
    where R is the nuclear radius and t a thickness parameter (default 2.3).
    Density is normalised to Z:
     \f[ \int rho(r) dr = Z \f]
    From this a potential is created that is normalised to Z/r as r->infinity.
 */
class NucleusDecorator: public LocalPotentialDecorator<NucleusDecorator>
{
public:
    NucleusDecorator(pHFOperator wrapped_hf, pOPIntegrator integration_strategy = pOPIntegrator());
    NucleusDecorator(const NucleusDecorator& other):
        LocalPotentialDecorator(other), nuclear_radius(other.nuclear_radius), nuclear_thickness(other.nuclear_thickness)
    {}

    /** Set parameters of proton distribution (in Fermi). */
    virtual void SetFermiParameters(double radius, double thickness = 2.3);
    virtual double GetNuclearRadius() const    { return nuclear_radius;    }    //!< Radius of nucleus in Fermi.
    virtual double GetNuclearThickness() const { return nuclear_thickness; }    //!< Nuclear thickness in Fermi.

    /** Get the nuclear density function, rho(r). */
    virtual RadialFunction GetNuclearDensity() const;

    /** Calculate nuclear RMS radius for the distribution used. */
    virtual double CalculateNuclearRMSRadius() const;

protected:
    virtual RadialFunction CalculateNuclearDensity(double radius, double thickness) const;

protected:
    double nuclear_radius;      //!< In Fermi
    double nuclear_thickness;   //!< In Fermi
};

typedef std::shared_ptr<NucleusDecorator> pNucleusDecorator;
typedef std::shared_ptr<const NucleusDecorator> pNucleusDecoratorConst;

#endif
