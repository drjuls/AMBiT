#ifndef SIGMA_POTENTIAL_H
#define SIGMA_POTENTIAL_H

#include "Universal/SpinorFunction.h"
#include "Universal/Lattice.h"
#include <vector>
#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> SigmaMatrix;

/** Storage class for Brueckner "Sigma" potential made of four radial matrices:
        Sigma(r1, r2) = ( ff  fg )
                        ( gf  gg )
    where each quadrant is a matrix in (r1, r2).
    By default it only stores one matrix: the "ff" part; use IncludeLower() to use the "fg" or "gg" parts.
 */
class SigmaPotential
{
public:
    SigmaPotential(pLattice lattice);
    /** Matrix size = (end_point - start_point)/stride */
    SigmaPotential(pLattice lattice, unsigned int end_point, unsigned int start_point = 0, unsigned int stride = 4);
    ~SigmaPotential();

    /** include_fg: store off diagonal (fg and gf) terms;
        include_gg: store gg term (can't imagine that it would be useful; smaller than ff by a factor ~(Z_ion * alpha)^2)
     */
    void IncludeLower(bool include_fg = false, bool include_gg = false);

    unsigned int size() const;
    void clear();
    void resize_and_clear(unsigned int new_size);

    /** Sigma(r1, r2) += s1(r1) * s2(r2) * coeff
        PRE: s1.size() & s2.size() >= size()
     */
    void AddToSigma(const SpinorFunction& s1, const SpinorFunction& s2, double coeff);

    /** Sigma.ff(r1, r2) += f1(r1) * f2(r2) * coeff
        PRE: s1.size() & s2.size() >= size()
     */
    void AddToSigma(const std::vector<double>& f1, const std::vector<double>& f2, double coeff);

    /** Return Integral[ Sigma(r1, r2). a(r2). dr2].
        PRE: a.size() >= size()
        Direct integration is hard-coded here for efficiency when using Eigen, therefore lattice is
        needed rather than an integrator.
     */
    SpinorFunction ApplyTo(const SpinorFunction& a) const;

    /** Attempt to read file. Return false if file not found, in which case SigmaPotential is not changed. */
    bool Read(const std::string& filename);
    void Write(const std::string& filename) const;

protected:
    // A matrix for each quadrant of Sigma
    pLattice lattice;
    SigmaMatrix ff, fg, gf, gg;
    bool use_fg, use_gg;

    unsigned int start;
    unsigned int matrix_size;
    unsigned int stride;
    std::vector<double> Rgrid;
    std::vector<double> dRgrid;
};

typedef std::shared_ptr<SigmaPotential> pSigmaPotential;

#endif
