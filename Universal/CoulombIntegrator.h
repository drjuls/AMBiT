#ifndef COULOMB_INTEGRATOR_H
#define COULOMB_INTEGRATOR_H

#include "Integrator.h"
#include "HartreeFock/State.h"

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
    CoulombIntegrator(Lattice& lat): Integrator(lat) {}
    virtual ~CoulombIntegrator(void) {}

    /** Calculate k-th term of coulomb potential due to charge density.
        if(charge != 0)
            Renormalise the density to charge.
     */
    void CoulombIntegrate(std::vector<double>& density, std::vector<double>& potential, unsigned int k, double charge = 0.);

    /** Same thing, but faster and less accurate.
        Good for just calculating a matrix element or something similar.
        There is no option to renormalise the density, because this wouldn't be very accurate anyway.
     */
    void FastCoulombIntegrate(const std::vector<double>& density, std::vector<double>& potential, unsigned int k);

    /** Get isotope shift between two states. In the first of these functions,
            f = s1.f
            L = s1.L
     */
    double IsotopeShiftIntegral(const std::vector<double> f, unsigned int L, const State& s2, std::vector<double>* P = NULL);
    double IsotopeShiftIntegral(const State& s1, const State& s2, std::vector<double>* P = NULL);

protected:
    class CoulombFunction : public Function2
    {
    public:
        CoulombFunction(Lattice& lat, const std::vector<double>& density, unsigned int k):
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
        Lattice& lattice;
        bool FirstPass;
    };
};

#endif