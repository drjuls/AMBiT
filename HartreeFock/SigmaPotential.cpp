#include "Include.h"
#include "SigmaPotential.h"

#define gap 4

SigmaPotential::SigmaPotential(pLattice lat, const std::string& file, unsigned int matrix_size, unsigned int start_point):
    lattice(lat), filename(file), changed_since_store(true), start(start_point)
{
    if(start%gap != 0)
        start = (unsigned int)(start/gap) * gap;

    size = (matrix_size - start)/gap;
    function = (double**)malloc(size * sizeof(double*));

    unsigned int i;
    for(i=0; i<size; i++)
    {   function[i] = (double*)malloc(size * sizeof(double));
        memset(function[i], 0, size * sizeof(double));
    }

    FILE* fp;
    if((fp = fopen(filename.c_str(), "rb")))
    {   fclose(fp);
        ReadIn();
    }
}

SigmaPotential::~SigmaPotential()
{
    for(unsigned int i=0; i<size; i++)
        free(function[i]);
    free(function);
}

std::vector<double> SigmaPotential::GetPotential(const std::vector<double>& f) const
{
    std::vector<double> ret(f.size());
    unsigned int limit = mmin((unsigned int)(f.size() - start)/gap, size);

    // get points from sigma function
    unsigned int i, j;
    const double* dR = lattice->dR();

    for(i = 0; i < limit; i++)
    {   // (j < i) swap indices (symmetry)
        for(j = 0; j < i; j++)
        {   ret[i*gap + start] += function[j][i] * f[j*gap + start] * dR[j*gap + start] * gap;
        }

        for(j = i; j < limit; j++)
        {   ret[i*gap + start] += function[i][j] * f[j*gap + start] * dR[j*gap + start] * gap;
        }
    }

    if(gap == 1)
        return ret;

    // otherwise gap == 4
    for(i = start; i<ret.size(); i+=gap)
    {   ret[i] = ret[i] * dR[i];
    }

    // interpolate for intermediate points
    limit = (limit - 1) * gap + start;

    // start points
    double A = ret[start];
    double B = ret[start + 4];
    double C = ret[start + 8];
    double D;
    ret[start + 1] = (21.*A + 14.*B - 3.*C)/32.;
    ret[start + 2] = (3.*A + 6.*B - C)/8.;
    ret[start + 3] = (5.*A + 30.*B - 3.*C)/32.;

    // end points
    A = ret[limit - 8];
    B = ret[limit - 4];
    C = ret[limit];
    ret[limit - 3] = (-3.*A + 30.*B + 5.*C)/32.;
    ret[limit - 2] = (-A + 6.*B + 3.*C)/8.;
    ret[limit - 1] = (-3.*A + 14.*B + 21.*C)/32.;

    // middle points
    for(j = start + 4; j <= limit - 8; j+=4)
    {
        A = ret[j - 4];
        B = ret[j];
        C = ret[j + 4];
        D = ret[j + 8];
        ret[j + 1] = (-7.*A + 105.*B + 35.*C - 5.*D)/128.;
        ret[j + 2] = (-A + 9.*B + 9.*C - D)/16.;
        ret[j + 3] = (-5.*A + 35.*B + 105.*C - 7.*D)/128.;
    }

    for(i=start; i<ret.size(); i++)
        ret[i] = ret[i]/dR[i];

    return ret;
}

double SigmaPotential::GetMatrixElement(const std::vector<double>& f1, const std::vector<double>& f2) const
{
    std::vector<double> pot = GetPotential(f2);

    const double* dR = lattice->dR();
    double value = 0.;
    for(unsigned int i=0; i < mmin(pot.size(), f1.size()); i++)
        value += pot[i] * f1[i] * dR[i];

    return value;
}

double SigmaPotential::GetSigma(unsigned int r1, unsigned int r2)
{
    if((r1 >= size*gap + start) || (r2 >= size*gap + start) || (r1 < start) || (r2 < start))
        return 0.;

    if(r1 >= r2)
        return function[(r2 - start)/gap][(r1 - start)/gap];
    else
        return function[(r1 - start)/gap][(r2 - start)/gap];
}

void SigmaPotential::AddToSigma(const std::vector<double>& f1, const std::vector<double>& f2, double coeff)
{
    if(((f1.size() - start)/gap > size) || ((f2.size() - start)/gap > size))
    {   
        ReSize((mmax(f1.size(), f2.size()) - start)/gap);
    }

    for(unsigned int i = 0; i<(f1.size() - start)/gap; i++)
        for(unsigned int j = i; j<(f2.size() - start)/gap; j++)
        {
            function[i][j] = function[i][j] + f1[i*gap + start] * f2[j*gap + start] * coeff;
        }

    changed_since_store = true;
}

void SigmaPotential::AddToSigma(const std::vector<double>& f1, const std::vector<double>& f2, double coeff, unsigned int f1_size, unsigned int f2_size)
{
    if(((f1_size - start)/gap > size) || ((f2_size - start)/gap > size))
    {   
        ReSize((mmax(f1_size, f2_size) - start)/gap);
    }

    for(unsigned int i = 0; i<(f1_size - start)/gap; i++)
        for(unsigned int j = i; j<(f2_size - start)/gap; j++)
        {
            function[i][j] = function[i][j] + f1[i*gap + start] * f2[j*gap + start] * coeff;
        }

    changed_since_store = true;
}

void SigmaPotential::AddDiagonal(const std::vector<double>& f, double coeff)
{
    unsigned int limit = (f.size() - start)/gap;
    limit = mmin(limit, size);

    for(unsigned int i = 0; i<limit; i++)
       function[i][i] = function[i][i] + f[i*gap + start] * coeff;

    changed_since_store = true;
}

void SigmaPotential::Store() const
{
    if(changed_since_store)
        WriteOut();
}

void SigmaPotential::ReadIn()
{
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        exit(1);

    unsigned int new_size;
    fread(&new_size, sizeof(unsigned int), 1, fp);
    fread(&start, sizeof(unsigned int), 1, fp);

    ReSize(new_size);

    for(unsigned int i = 0; i<size; i++)
    {
        fread(function[i], sizeof(double), size, fp);
    }

    fclose(fp);

    changed_since_store = false;
}

void SigmaPotential::WriteOut() const
{
    FILE* fp = fopen(filename.c_str(), "wb");

    if(!fp)
        exit(1);

    fwrite(&size, sizeof(unsigned int), 1, fp);
    fwrite(&start, sizeof(unsigned int), 1, fp);

    for(unsigned int i = 0; i<size; i++)
    {
        fwrite(function[i], sizeof(double), size, fp);
    }

    fclose(fp);

    changed_since_store = false;
}

unsigned int SigmaPotential::Size()
{
    return (size * gap + start);
}

void SigmaPotential::ReSize(unsigned int new_size)
{
    if(new_size < size)
    {
        *logstream << "  Sigma resize: " << size << " -> " << new_size << std::endl;
    	unsigned int i;
	    for(i = new_size; i<size; i++)
	        free(function[i]);

	    function = (double**)realloc(function, new_size * sizeof(double*));

        for(i = 0; i<new_size; i++)
	        function[i] = (double*)realloc(function[i], new_size * sizeof(double));
    }

    else if(new_size > size)
    {
        *logstream << "  Sigma resize: " << size << " -> " << new_size << std::endl;
        function = (double**)realloc(function, new_size * sizeof(double*));

        unsigned int i, j;
        for(i=0; i<size; i++)
            function[i] = (double*)realloc(function[i], new_size * sizeof(double));

        for(i=size; i<new_size; i++)
        {   function[i] = (double*)malloc(new_size * sizeof(double));
            memset(function[i], 0, new_size * sizeof(double));
        }        

        for(i=0; i<size; i++)
	    {   for(j=size; j<new_size; j++)
	            function[i][j] = 0.;
	    }
    }

    size = new_size;
}

void SigmaPotential::Reset()
{
  for(unsigned int i=0; i<size; i++)
      memset(function[i], 0, size * sizeof(double));
}
