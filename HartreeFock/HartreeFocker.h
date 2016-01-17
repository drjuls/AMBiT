#ifndef HARTREE_FOCKER_H
#define HARTREE_FOCKER_H

#include "HFOperator.h"
#include "Universal/Enums.h"

class HartreeFocker
{
public:
    typedef std::function<double(HartreeFocker&, pOrbital, pHFOperator, pSpinorFunctionConst)> IteratorFunction;

public:
    HartreeFocker(pODESolver ode_solver): odesolver(ode_solver), continuum_normalisation_type(ContinuumNormalisation::LandauEnergy) {}

    /** Get first guess at the core orbitals using Thomas-Fermi and gradually building in
        a self-consistent direct potential.
        PRE: core should have all orbitals in place, albeit without radial functions.
     */
    void StartCore(pCore core, pHFOperator hf);

    /** Iterate all orbitals in core until self-consistency is reached. */
    void SolveCore(pCore core, pHFOperator hf);

    /** Create a new orbital in the field of the core.
        Return number of loops required for HF convergence.
     */
    unsigned int CalculateExcitedState(pOrbital orbital, pHFOperator hf);

    /** Create a new continuum wavefunction in the field of the core.
        Return number of loops required for HF convergence.
        Return value of 0 means start_sine was not reached and method failed.
     */
    unsigned int CalculateContinuumWave(pContinuumWave s, pHFOperator hf);

    /** Find energy eigenvalue for orbital with a given exchange potential.
        Note: this function does not iterate/update the exchange potential,
            so the final orbital is not an eigenvalue of the hf operator.
        Return change in energy.
     */
    double ConvergeOrbital(pOrbital orbital, pHFOperator hf, pSpinorFunctionConst exchange, IteratorFunction iterator = &HartreeFocker::IterateOrbitalTailMatching, double energy_tolerance = 1.e-8);

    /** Find self-consistent solution to hf operator, including exchange.
        Exchange is calculated from hf operator in each loop.
        Return change in energy.
     */
    double ConvergeOrbitalAndExchange(pOrbital orbital, pHFOperator hf, IteratorFunction iterator = &HartreeFocker::IterateOrbital, double energy_tolerance = 1.e-14);

    /** Iterator function for ConvergeOrbital and ConvergeOrbitalAndExchange.
        IterateOrbital uses the Greens method discussed in Johnson. A non-zero exchange term is required.
        Return change in energy.
     */
    double IterateOrbital(pOrbital orbital, pHFOperator hf, pSpinorFunctionConst exchange);

    /** Iterator function for ConvergeOrbital and ConvergeOrbitalAndExchange using tail matching method.
        Unlike the Greens method used in IterateOrbital(), this method doesn't require a
        source term (exchange term), so it can solve without exchange or with a local exchange approximation.
        Return change in energy.
     */
    double IterateOrbitalTailMatching(pOrbital orbital, pHFOperator hf, pSpinorFunctionConst exchange = pSpinorFunction());

    unsigned int IntegrateContinuum(pContinuumWave s, pHFOperator hf, pSpinorFunction exchange, double& final_amplitude, double& final_phase);

    double WavefunctionTolerance = 1.e-11;
    double EnergyTolerance = 1.e-14;
    double TailMatchingEnergyTolerance = 1.e-8;
    ContinuumNormalisation continuum_normalisation_type;

protected:
    pODESolver odesolver;
    unsigned int MaxHFIterations = 500;

protected:
    class ContinuumPhaseODE: public OneDimensionalODE
    {
    public:
        ContinuumPhaseODE(pLattice lattice, const RadialFunction& hf_direct, int kappa, double energy):
            OneDimensionalODE(lattice), V(hf_direct), kappa(kappa), energy(energy)
        {}

        /** Get df/dr = w[0] given point r, f.
            PRE: w should be an allocated double.
         */
        virtual void GetODEFunction(unsigned int latticepoint, const RadialFunction& f, double* w) const
        {
            double R = lattice->R(latticepoint);
            *w = sqrt(2. * fabs(energy + V.f[latticepoint] - double(kappa*(kappa+1))/(2. * R * R)));
        }

        /** Get numerical coefficients of the ODE at the point r, f.
            PRE: w_f and w_const should be allocated 2 dimensional arrays.
         */
        virtual void GetODECoefficients(unsigned int latticepoint, const RadialFunction& f, double* w_f, double* w_const) const
        {
            *w_f = 0.;
            GetODEFunction(latticepoint, f, w_const);
        }

        /** Get Jacobian dw[i]/df and dw[i]/dr at a point r, f.
            PRE: jacobian and dwdr should allocated doubles.
         */
        virtual void GetODEJacobian(unsigned int latticepoint, const RadialFunction& f, double* jacobian, double* dwdr) const
        {
            *jacobian = 0.;

            double R = lattice->R(latticepoint);
            double w;
            GetODEFunction(latticepoint, f, &w);
            *dwdr = 1./w * (V.dfdr[latticepoint] + double(kappa*(kappa+1))/(R*R*R));
        }

        /** Get approximation to solution for first numpoints near the origin. */
        virtual void EstimateSolutionNearOrigin(unsigned int numpoints, RadialFunction& f) const
        {
            double w_i;
            const double* dR = lattice->dR();

            unsigned int i = 0;

            GetODEFunction(i, f, &w_i);
            f.f[i] = 0.5 * w_i * dR[i];
            f.dfdr[i] = w_i;

            for(i=1; i < numpoints; i++)
            {
                GetODEFunction(i, f, &w_i);
                f.dfdr[i] = w_i;
                
                f.f[i] = f.f[i-1] + 0.5 * (f.dfdr[i-1] * dR[i-1] + f.dfdr[i] * dR[i]);
            }
        }

        /** Get approximation to solution for last numpoints far from the origin. */
        virtual void EstimateSolutionNearInfinity(unsigned int numpoints, RadialFunction& f) const
        {   exit(2);
        }
        
    protected:
        const RadialFunction& V;
        int kappa;
        double energy;
    };
};

#endif
