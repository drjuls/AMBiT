#include "Include.h"
#include "Lattice.h"
#include <fstream>

// Take a little endian, 4 byte DWORD and contract to a number
inline unsigned long DWORD_LE_CONTRACT(const char* dw)
{
    unsigned long ret =
        (static_cast<unsigned char>(dw[0])) +
        (static_cast<unsigned char>(dw[1])<<8) +
        (static_cast<unsigned char>(dw[2])<<16) +
        (static_cast<unsigned char>(dw[3])<<24);
    return ret;
}

void CollectFloatingValues(const char* buff, double* values, unsigned int size)
{
    unsigned long lvalue;
    float& fvalue = reinterpret_cast<float&>(lvalue);
    const char* pbuff = buff;

    unsigned int i = 0;
    while(i < size)
    {   lvalue = DWORD_LE_CONTRACT(pbuff);
        values[i] = fvalue;
        i++;
        pbuff += 4;
    }
}

Lattice::Lattice(const Lattice& other):
    beta(other.beta), NumPoints(other.NumPoints), rmin(other.rmin), h(other.h)
{
    r = (double*)malloc(NumPoints * sizeof(double));
    dr = (double*)malloc(NumPoints * sizeof(double));
    
    for(unsigned int i=0; i<NumPoints; i++)
    {   r[i] = other.r[i];
        dr[i] = other.dr[i];
    }

}

Lattice::Lattice(unsigned int numpoints, double r_min, double r_max):
    beta(4.0), NumPoints(numpoints), rmin(r_min)
{
    r = (double*)malloc(NumPoints * sizeof(double));
    dr = (double*)malloc(NumPoints * sizeof(double));
    
    h =(r_max - rmin + beta*log(r_max/rmin))/(NumPoints-1);
    
    for(unsigned int i=0; i<NumPoints; i++)
    {   r[i] = lattice_to_real(i);
        dr[i] = calculate_dr(r[i]);
    }
}

Lattice::Lattice(const std::string& filename)
{
    std::fstream infile;
    infile.open(filename.c_str(), std::ios::in | std::ios::binary);
    
    if(!infile)
        return;

    unsigned int size = 1000;
    r = new double[size];
    dr = new double[size];
    NumPoints = size;

    char short_buff[4];
    infile.clear();
    infile.read(short_buff, 4);
    unsigned int block_size = static_cast<unsigned int>(DWORD_LE_CONTRACT(short_buff));

    char* buffer = new char[block_size];

    // Skip potential
    infile.read(buffer, 4*size);

    // Get lattice definition
    infile.read(buffer, 4*size);
    CollectFloatingValues(buffer, r, size);

    infile.read(buffer, 4*size);
    CollectFloatingValues(buffer, dr, size);

    infile.read(buffer, 4*181);
    infile.read(short_buff, 4);

    beta = 4;

    unsigned long lvalue = DWORD_LE_CONTRACT(short_buff);
    float& fvalue = reinterpret_cast<float&>(lvalue);
    h = fvalue;

    infile.close();

    for(unsigned int i=0; i<size; i++)
    {   dr[i] = (dr[i]) * (r[i]) * h;
    }
    rmin = r[0];

    delete[] buffer;
}

Lattice::~Lattice(void)
{   free(r);
    free(dr);

    std::vector<double*>::iterator it = r_power.begin();
    while(it != r_power.end())
    {   free(*it);
        it++;
    }
}

void Lattice::ReSize(double min_size)
{
    if(min_size <= r[NumPoints-1])
        return;

    unsigned int old_size = NumPoints;

    do
    {   NumPoints *= 2;
    }while(lattice_to_real(NumPoints - 1) <= min_size);

    r = (double*)realloc(r, NumPoints * sizeof(double));
    dr = (double*)realloc(dr, NumPoints * sizeof(double));

    for(unsigned int i=old_size; i<NumPoints; i++)
    {   r[i] = lattice_to_real(i);
        dr[i] = calculate_dr(r[i]);
    }

    // Do powers of r
    for(unsigned int k =0; k < r_power.size(); k++)
    {
        r_power[k] = (double*)realloc(r_power[k], NumPoints * sizeof(double));
        double* points = r_power[k];

        memcpy(points+old_size, r+old_size, (NumPoints - old_size) * sizeof(double));
        
        for(unsigned int i=old_size; i<NumPoints; i++)
        {
            for(unsigned int j = 2; j<k+2; j++)
                points[i] = points[i] * r[i];
        }
    }
}

double Lattice::lattice_to_real(unsigned int i) const
{
    // Lattice coordinate x = rmin + ih,
    //      y = x + beta * log(rmin)
    double y = rmin + beta*log(rmin) + h*double(i);
    double r, rold;

    if(beta <= y)
        r = y;
    else
        r = exp(y/beta);

    // Solve using Newton's method f(r) == 0, where
    //    f(r) = r + beta*log(r/rmin) - x
    //         = r + beta*log(r) - y
    do
    {   rold = r;
        r = r - (r + beta*log(r) - y)/(1. + beta/r);
        if(r <= 0.0)
            r = rold/2.0;
    } while(fabs(r - rold)/r > 1.e-13);

    return r;
}

unsigned int Lattice::real_to_lattice(double r_point)
{
    if(r_point > MaxRealDistance())
        ReSize(r_point);
    else if (r_point <= r[0])
        return 0;

    unsigned int i_min = 0;
    unsigned int i_max = NumPoints;
    unsigned int i_mid;

    // Bisection search
    while(i_max - i_min > 1)
    {   i_mid = (i_min + i_max)/2;
        (r_point >= r[i_mid]) ? i_min = i_mid : i_max = i_mid;
    }

    if(r_point == r[i_min])
        return i_min;
    else
        return i_min + 1;
}

bool Lattice::operator==(const Lattice& other) const
{
    return ((beta == other.beta) && (h == other.h) && (rmin == other.rmin));
}

double Lattice::calculate_dr(double r_point) const
{
    return r_point/(beta + r_point) * h;
}

const double* Lattice::Calculate_Rpower(unsigned int k)
{
    // Create all powers up to and including k
    for(unsigned int power = r_power.size()+2; power <= k; power++)
    {
        double* points = (double*)malloc(NumPoints * sizeof(double));
        if(power == 2)
            memcpy(points, r, NumPoints * sizeof(double));
        else
            memcpy(points, r_power[power-3], NumPoints * sizeof(double));

        for(unsigned int i=0; i<NumPoints; i++)
        {
            points[i] = points[i] * r[i];
        }

        r_power.push_back(points);
    }

    return r_power[k-2];
}
