#ifndef SIGMA_POTENTIAL_H
#define SIGMA_POTENTIAL_H

#include "Universal/Lattice.h"
#include <vector>

class SigmaPotential
{
public:
    SigmaPotential(Lattice* lat, const std::string& file, unsigned int matrix_size = 0, unsigned int start_point = 0);
    ~SigmaPotential();

    std::vector<double> GetPotential(const std::vector<double>& f) const;

    double GetSigma(unsigned int r1, unsigned int r2);

    void AddToSigma(const std::vector<double>& f1, const std::vector<double>& f2, double coeff);

    unsigned int Size() { return size; }
    void ReSize(unsigned int new_size);

    void Reset();

    void Store() const;

protected:
    void WriteOut() const;
    void ReadIn();

    mutable bool changed_since_store;
    std::string filename;

    unsigned int size;
    unsigned int start;
    double** function;

    Lattice* lattice;
};

#endif
