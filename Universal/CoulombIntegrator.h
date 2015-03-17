#ifndef COULOMB_INTEGRATOR_H
#define COULOMB_INTEGRATOR_H

#include "Integrator.h"
#include "HartreeFock/Orbital.h"

class CoulombIntegrator : public Integrator
{
    /** Find the potential due to a charge density
                              r< ^k
            I(r) = Integral[ ---------.density(r').dr' ]
                             r> ^(k+1)
        Split into two parts,
                      1
           I1(r) = ------- . Integral[ r' ^k .density(r').dr' ] 
                   r ^(k+1)    0->r
                                          1
            I2(r) = r ^k . Integral  [ ------- .density(r').dr' ] 
                          r->infinity  r'^(k+1)

          dI1/dr = -(k+1)/r .I1 + density/r
          dI2/dr =    k/r .I2   - density/r
     */
public:
    CoulombIntegrator(pLattice lat): Integrator(lat) {}
    virtual ~CoulombIntegrator(void) {}

    /** Calculate k-th term of coulomb potential due to charge density.
        if(charge != 0)
            Renormalise the density to charge.
     */
    void CoulombIntegrate(std::vector<double>& density, std::vector<double>& potential, unsigned int k, double charge = 0.);

    /** Same thing, but faster and less accurate.
        Good for just calculating a matrix element or something similar.
        There is no option to renormalise the density, because this wouldn't be very accurate anyway.
        POST: potential.size = density.size
     */
    void FastCoulombIntegrate(const std::vector<double>& density, std::vector<double>& potential, unsigned int k);

    /** This one doesn't change the potential size.
        PRE: density.size >= density_size (after which all density points are assumed zero).
             lattice.size >= potential_size
             Also, for accuracy, potential.size >= density_size.
     */
    void FastCoulombIntegrate(const std::vector<double>& density, std::vector<double>& potential, unsigned int k, unsigned int density_size);

protected:
    class CoulombFunction : public Function2
    {
    public:
        CoulombFunction(pLattice lat, const std::vector<double>& density, unsigned int k):
            Function2(), Density(&density), lattice(lat), K(k), FirstPass(true) {}
        virtual ~CoulombFunction(void) {}

        virtual double Coeff1(int point) const;
        virtual double Coeff2(int point) const;

        virtual std::vector<double> Coeff1(void) const;
        virtual std::vector<double> Coeff2(void) const;

        void SetDensity(const std::vector<double>& density) { Density = &density; }
        void SetDirection(bool first_pass) { FirstPass = first_pass; }
    protected:
        double K;
        const std::vector<double>* Density;
        pLattice lattice;
        bool FirstPass;
    };
};

#endif
