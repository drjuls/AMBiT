#ifndef CONVOLVE_DR_H
#define CONVOLVE_DR_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <new>
#include <sstream>
#include <map>

class Convolver
{
public:
    Convolver(): num_bins(0), bin_cross_sections(NULL) {}
    ~Convolver()
    {   if(bin_cross_sections)
            delete[] bin_cross_sections;
    }
  
    void ReadCrossSections(const std::string filename = "dr.txt");
    void BinCrossSections(double energy_width, double energy_max);
    void Convolve(double temp_perp, double temp_parallel);
    double ErrorFunction(double x);
    
protected:
    /** Map energy to integrated (partial) DR cross-section, as input from file. */
    std::map<double, double> IntegratedDRCrossSections;

    unsigned int num_bins;
    double bin_width;
    double* bin_cross_sections;
};

#endif