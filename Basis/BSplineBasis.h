#ifndef BSPLINE_BASIS_H
#define BSPLINE_BASIS_H

enum class SplineType {NotreDame, Reno, Vanderbilt};

typedef Orbital BSpline;
typedef pOrbital pBSpline;

/** Make basis by diagonalising hf operator over a set of B-splines generated using GSL.
    Can generate negative energy states (Dirac sea) using GenerateCompleteBasis().
    Spline types Reno and Vanderbilt implement dual kinetic balance;
    see Beloy and Derevianko, Comp. Phys. Comm. 179, 310 (2008).
 */
class BSplineBasis
{
public:
    BSplineBasis(pLattice lattice, int n, int k, double rmax, double dr0 = 0.0, SplineType method = SplineType::Reno):
        lattice(lattice)
    {   SetParameters(n, k, rmax, dr0, method);
    }
    virtual ~BSplineBasis() {}

    void SetParameters(int n, int k, double rmax, SplineType method = SplineType::Reno);
    void SetParameters(int n, int k, double rmax, double dr0, SplineType method = SplineType::Reno);

    /** Generate B-splines from core states up to max_pqn. */
    pOrbitalMap GenerateBSplines(pHFOperator hf, int kappa, int max_pqn);

    /** Return all positive energy states (removes Dirac sea). */
    pOrbitalMap GeneratePositiveBasis(pHFOperator hf, int kappa);

    /** Generate complete B-spline basis set including negative energy states (Dirac sea). */
    pOrbitalMap GenerateCompleteBasis(pHFOperator hf, int kappa);

protected:
    /** Generate a set of splines using GSL routines. */
    std::vector<pBSpline> MakeSplines(int kappa, pPhysicalConstant constants);

protected:
    pLattice lattice;
    int N, K;
    double Rmax, dR0;
    SplineType spline_type;
};

#endif
