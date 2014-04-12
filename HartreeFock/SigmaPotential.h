#ifndef SIGMA_POTENTIAL_H
#define SIGMA_POTENTIAL_H

#include "Universal/Lattice.h"
#include <vector>

class SigmaPotential
{
public:
    SigmaPotential(pLattice lat, const std::string& file = "", unsigned int matrix_size = 0, unsigned int start_point = 0);
    ~SigmaPotential();

    /** Get potential
            V(r1) = Integral[Sigma(r1, r2) * f(r2), {r2, 0, Infinity}]
     */
    std::vector<double> GetPotential(const std::vector<double>& f) const;

    /** Get matrix element <f1 | Sigma | f2> */
    double GetMatrixElement(const std::vector<double>& f1, const std::vector<double>& f2) const;

    double GetSigma(unsigned int r1, unsigned int r2);

    /** function(r1, r2) += f1(r1) * f2(r2) * coeff
        POST: Sigma.size() >= f1.size() and f2.size()
     */
    void AddToSigma(const std::vector<double>& f1, const std::vector<double>& f2, double coeff);

    /** PRE:  f1.size() >= f1_size, f2.size() >= f2_size
        POST: Sigma.size() >= f1_size and f2_size
     */
    void AddToSigma(const std::vector<double>& f1, const std::vector<double>& f2, double coeff, unsigned int f1_size, unsigned int f2_size);

    /** function(r1, r1) += f(r1) * coeff */
    void AddDiagonal(const std::vector<double>& f, double coeff);
    
    unsigned int Size();

    void Reset();

    void Store() const;

protected:
    void WriteOut() const;
    void ReadIn();

    void Size(unsigned int new_size);
    
    mutable bool changed_since_store;
    std::string filename;

    unsigned int size;  // Physical size of function array
    unsigned int start;
    double** function;

    pLattice lattice;
};

#endif
