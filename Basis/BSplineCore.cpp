#include "BSplineCore.h"
#include "Include.h"
#include "Universal/Interpolator.h"
#include "Universal/CoulombIntegrator.h"
#include "Spline.h"

#define NON_REL_SCALING true

BSplineCore::BSplineCore(BSplineGrid* lat, const Core* core):
    hfcore(core), grid(lat)
{
    order = 6;
    Lattice* hflattice = hfcore->GetLattice();

    const double* x = grid->R();
    const double* dx = grid->dR();

    // Interpolate potential
    // Multiply by r to smooth potential before interpolation
    // TODO: Remove multiplication by r at small r (V is almost constant here)?
    std::vector<double> hfpot = hfcore->GetHFPotential();
    const double* R = hflattice->R();
    for(unsigned int i=0; i<hfpot.size(); i++)
        hfpot[i] = hfpot[i] * R[i];

    potential.clear();
    potential.resize(grid->Size());
    unsigned int point;

    Interpolator interpHFLat(hflattice);
    double V, dV;
    for(point = 0; point < grid->Size(); point++)
    {   
        interpHFLat.Interpolate(hfpot, x[point], V, dV, order);
        // divide by x again
        potential[point] = V/x[point];
    }
}
/*
FROM CONSTRUCTOR:

    // Interpolate wavefunctions
    const DiscreteState* hfstate;
    DiscreteState* s;
    
    ConstStateIterator it = hfcore->GetConstStateIterator();
    while(!it.AtEnd())
    {
        hfstate = it.GetState();
        s = new DiscreteState(grid, hfstate->RequiredPQN(), hfstate->Kappa());
        s->SetOccupancy(hfstate->Occupancy());
        s->SetNu(hfstate->Nu());

        double realsize = hflattice->R(hfstate->Size()-1);
        unsigned int size = grid->real_to_lattice(realsize);
        s->ReSize(size);

        double dfdx, dgdx;
        for(point = 0; point < size; point++)
        {   
            interpHFLat.Interpolate(hfstate->f, x[point], s->f[point], dfdx, order);
            interpHFLat.Interpolate(hfstate->g, x[point], s->g[point], dgdx, order);
            s->df[point] = dfdx * dx[point];
            s->dg[point] = dgdx * dx[point];
        }

        AddState(s);
        it.Next();
    }

void BSplineCore::CalculateExchange(unsigned int bspline, int kappa, bool upper, CoupledFunction& exchange)
{
    unsigned int L;
    if(kappa > 0)
        L = (unsigned int)(kappa);
    else
        L = (unsigned int)(-kappa-1);

    double J = (double)abs(kappa) - 0.5;

    // Lattices
    Lattice* hflattice = hfcore->GetLattice();
    const double* R = hflattice->R();
    Interpolator spline_interp(hflattice);
    CoulombIntegrator HFI(*hflattice);
    //CoupledFunction hfexchange(hflattice->Size());

    const double* X = grid->R();
    exchange.Clear();
    exchange.ReSize(grid->Size());

    // Get BSpline on HFlattice
    std::vector<double> hfspline(hflattice->Size());

    const double* knots = grid->GetSplineKnots();
    int k = grid->GetK();
    int nderiv = 1;
    double* fspline_buf = new double[nderiv * MaximumK];
    int left = (int)bspline;

    unsigned int b_lower_point = hflattice->real_to_lattice(knots[bspline]);    // check!!!
    unsigned int b_upper_point = hflattice->real_to_lattice(knots[bspline+k]);
    unsigned int i;

    for(i = b_lower_point; i < b_upper_point; i++)
    {   
        double r = R[i];
        while(r > knots[left+1])
            left++;

        int leftplus1 = left + 1;   // Fortran array index starts from 1
        bsplvd(knots, &k, &r, &leftplus1, fspline_buf, &nderiv);
        hfspline[i] = fspline_buf[bspline + k - left - 1];
    }

    delete[] fspline_buf;

    // Get BSpline on spline lattice
    //std::vector<double> xspline(grid->Size());
    //unsigned int spline_lower_point = grid->real_to_lattice(knots[bspline]);
    //unsigned int spline_upper_point = grid->real_to_lattice(knots[bspline+k]);

    //for(i = spline_lower_point; i < spline_upper_point; i++)
    //{
    //    double x = X[i];
    //    while(x > knots[left+1])
    //        left++;

    //    int leftplus1 = left + 1;   // Fortran array index starts from 1
    //    bsplvd(knots, &k, &x, &leftplus1, fspline_buf, &nderiv);
    //    xspline[i] = fspline_buf[bspline + k - left - 1];
    //}

    // Sum over all core states
    ConstDiscreteStateIterator cs = hfcore->GetConstDiscreteStateIterator();
    while(!cs.AtEnd())
    {
        const DiscreteState& other = cs.GetState();
        const State* xother = GetState(StateInfo(&other));

        unsigned int lower_point = b_lower_point;
        unsigned int upper_point = other.Size();

        // Get overlap of wavefunctions
        std::vector<double> density(hflattice->Size());
        if(upper)
            for(unsigned int i=lower_point; i < upper_point; i++)
                density[i] = other.f[i] * hfspline[i];
        else
            for(unsigned int i=lower_point; i < upper_point; i++)
                density[i] = Constant::AlphaSquared * other.g[i] * hfspline[i];

        // Sum over all k
        for(unsigned int k = abs((int)other.L() - (int)L); k <= (other.L() + L); k+=2)
        {
            double coefficient = Constant::Wigner3j(k, J, other.J(), 0., .5, -.5);
            coefficient = (2 * abs(other.Kappa())) * coefficient * coefficient;

            // Open shells need to be scaled
            if(other.Occupancy() != double(2 * abs(other.Kappa())))
            {
                double ex = 1.;
                if(NON_REL_SCALING)
                {   // Average over non-relativistic configurations
                    if(other.Kappa() == -1)
                    {
                        if(kappa != other.Kappa())
                            ex = other.Occupancy()/double(2 * abs(other.Kappa()));
                        else if(k)
                            ex = (other.Occupancy()-1.)/double(2 * abs(other.Kappa()) - 1);
                    }
                    else
                    {
                        int other_kappa = - other.Kappa() - 1;
                        const DiscreteState* ds = hfcore->GetState(StateInfo(other.PQN(), other_kappa));

                        if((kappa != other.Kappa()) && (kappa != ds->Kappa()))
                            ex = (other.Occupancy() + ds->Occupancy())/double(2 * (abs(other.Kappa()) + abs(ds->Kappa())));
                        else if(k)
                            ex = (other.Occupancy() + ds->Occupancy() - 1.)/double(2 * (abs(other.Kappa()) + abs(ds->Kappa())) - 1);
                    }
                }
                else
                {   // Average over relativistic configurations
                    if(kappa != other.Kappa())
                        ex = other.Occupancy()/double(2 * (abs(other.Kappa())));
                    else if(k)
                        ex = (other.Occupancy() - 1.)/double(2 * (abs(other.Kappa())) - 1);
                }

                coefficient = coefficient * ex;
            }

            // Integrate density to get (1/r)Y(ab,r)
            std::vector<double> hfpotential;
            HFI.CoulombIntegrate(density, hfpotential, k);

            const double* hfRk = hflattice->Rpower(k);  // NULL if (k == 0)

            // Multiply potential by r^k and interpolate.
            if(hfRk)
            {   for(i=0; i<hfpotential.size(); i++)
                    hfpotential[i] = hfpotential[i] * hfRk[i];
            }

            std::vector<double> xpotential(xother->Size());
            double dpot;
            for(i=0; i<xpotential.size(); i++)
                spline_interp.Interpolate(hfpotential, X[i], xpotential[i], dpot, order);
            
            // Then create exchange in spline grid.
            const double* spXk = grid->Rpower(k);
            if(spXk)
            {   for(i=0; i<xpotential.size(); i++)
                    xpotential[i] = xpotential[i]/spXk[i];
            }

            for(i=0; i<xother->Size(); i++)
            {   exchange.f[i] = exchange.f[i] + coefficient * xpotential[i] * xother->f[i];
                exchange.g[i] = exchange.g[i] + coefficient * xpotential[i] * xother->g[i];
            }

            //if(hfcore->GetNuclearInverseMass() && (k == 1))
            //{
            //    std::vector<double> P(upper_point);
            //    double sms = HFI.IsotopeShiftIntegral(current, other, &P);
            //    
            //    for(unsigned int i=0; i<upper_point; i++)
            //    {
            //        hfexchange.f[i] = hfexchange.f[i] + coefficient * hfcore->GetNuclearInverseMass() * sms *  P[i];
            //    }
            //}
        }
        cs.Next();
    }
}
*/
void BSplineCore::CalculateExchange(unsigned int bspline, int kappa, bool upper, CoupledFunction& exchange)
{
    unsigned int L;
    if(kappa > 0)
        L = (unsigned int)(kappa);
    else
        L = (unsigned int)(-kappa-1);

    double J = (double)abs(kappa) - 0.5;

    // Lattices
    Lattice* hflattice = hfcore->GetLattice();
    const double* R = hflattice->R();

    CoulombIntegrator HFI(*hflattice);

    CoupledFunction hfexchange(hflattice->Size());

    // Get BSpline on HFlattice
    std::vector<double> spline(hflattice->Size());

    const double* knots = grid->GetSplineKnots();
    int k = grid->GetK();
    int nderiv = 1;
    double* fspline_buf = new double[nderiv * MaximumK];
    int left = (int)bspline;

    unsigned int b_lower_point = hflattice->real_to_lattice(knots[bspline]);    // check!!!
    unsigned int b_upper_point = hflattice->real_to_lattice(knots[bspline+k]);
    unsigned int i;

    for(i = b_lower_point; i < b_upper_point; i++)
    {   
        double r = R[i];
        while(r > knots[left+1])
            left++;

        int leftplus1 = left + 1;   // Fortran array index starts from 1
        bsplvd_(knots, &k, &r, &leftplus1, fspline_buf, &nderiv);
        spline[i] = fspline_buf[bspline + k - left - 1];
    }

    delete[] fspline_buf;

    // Sum over all core states
    ConstStateIterator cs = hfcore->GetConstStateIterator();
    while(!cs.AtEnd())
    {
        const DiscreteState& other = *cs.GetState();
        unsigned int lower_point = b_lower_point;
        unsigned int upper_point = other.Size();

        // Get overlap of wavefunctions
        std::vector<double> density(upper_point);
        if(upper)
            for(unsigned int i=lower_point; i < upper_point; i++)
                density[i] = other.f[i] * spline[i];
        else
            for(unsigned int i=lower_point; i < upper_point; i++)
                density[i] = other.g[i] * spline[i];

        // Sum over all k
        for(unsigned int k = abs((int)other.L() - (int)L); k <= (other.L() + L); k+=2)
        {
            double coefficient = Constant::Wigner3j(k, J, other.J(), 0., .5, -.5);
            coefficient = (2 * abs(other.Kappa())) * coefficient * coefficient;

            // Open shells need to be scaled
            if(other.Occupancy() != double(2 * abs(other.Kappa())))
            {
                double ex = 1.;
                if(NON_REL_SCALING)
                {   // Average over non-relativistic configurations
                    if(other.Kappa() == -1)
                    {
                        if(kappa != other.Kappa())
                            ex = other.Occupancy()/double(2 * abs(other.Kappa()));
                        else if(k)
                            ex = (other.Occupancy()-1.)/double(2 * abs(other.Kappa()) - 1);
                    }
                    else
                    {
                        int other_kappa = - other.Kappa() - 1;
                        const DiscreteState* ds = hfcore->GetState(StateInfo(other.RequiredPQN(), other_kappa));

                        if((kappa != other.Kappa()) && (kappa != ds->Kappa()))
                            ex = (other.Occupancy() + ds->Occupancy())/double(2 * (abs(other.Kappa()) + abs(ds->Kappa())));
                        else if(k)
                            ex = (other.Occupancy() + ds->Occupancy() - 1.)/double(2 * (abs(other.Kappa()) + abs(ds->Kappa())) - 1);
                    }
                }
                else
                {   // Average over relativistic configurations
                    if(kappa != other.Kappa())
                        ex = other.Occupancy()/double(2 * (abs(other.Kappa())));
                    else if(k)
                        ex = (other.Occupancy() - 1.)/double(2 * (abs(other.Kappa())) - 1);
                }

                coefficient = coefficient * ex;
            }

            // Integrate density to get (1/r)Y(ab,r)
            std::vector<double> potential;
            HFI.CoulombIntegrate(density, potential, k);

            // todo: Multiply potential by r^k and interpolate.
            //       Then create exchange in spline grid.

            for(unsigned int i=0; i<upper_point; i++)
            {
                hfexchange.f[i] = hfexchange.f[i] + coefficient * potential[i] * other.f[i];
                hfexchange.g[i] = hfexchange.g[i] + coefficient * potential[i] * other.g[i];
            }

            if(hfcore->GetNuclearInverseMass() && (k == 1))
            {
                std::vector<double> P(upper_point);
                double sms = HFI.IsotopeShiftIntegral(spline, L, other, &P);
                
                for(unsigned int i=0; i<upper_point; i++)
                {
                    hfexchange.f[i] = hfexchange.f[i] + coefficient * hfcore->GetNuclearInverseMass() * sms *  P[i];
                }
            }
        }
        cs.Next();
    }
    
    // Interpolate hfexchange onto spline grid
    exchange.Clear();
    exchange.ReSize(grid->Size());
    Interpolator spline_interp(hflattice);

    const double* x = grid->R();

    double deriv;
    for(i=0; i<grid->Size(); i++)
    {
        spline_interp.Interpolate(hfexchange.f, x[i], exchange.f[i], deriv, order);
        spline_interp.Interpolate(hfexchange.g, x[i], exchange.g[i], deriv, order);
    }
}
