#include "Include.h"
#include "SigmaPotential.h"
#include "Universal/Interpolator.h"

#define stride 4

#if (stride == 1)
    typedef Eigen::Map<const Eigen::VectorXd> EigenVectorMapped;
    typedef Eigen::Map<const Eigen::ArrayXd> EigenArrayMapped;
#else
    typedef Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<stride>> EigenVectorMapped;
    typedef Eigen::Map<const Eigen::ArrayXd, Eigen::Unaligned, Eigen::InnerStride<stride>> EigenArrayMapped;
#endif

SigmaPotential::SigmaPotential(pLattice lattice):
    start(0), matrix_size(0), lattice(lattice)
{}

SigmaPotential::SigmaPotential(pLattice lattice, unsigned int end_point, unsigned int start_point):
    start(start_point), matrix_size(0), lattice(lattice)
{
    resize_and_clear(end_point);
}

SigmaPotential::~SigmaPotential()
{}

void SigmaPotential::clear()
{
    ff.setZero();
    fg.setZero();
    gf.setZero();
    gg.setZero();
}

void SigmaPotential::resize_and_clear(unsigned int new_size)
{
    if(new_size > lattice->size())
        new_size = lattice->size();

    matrix_size = (new_size - start)/stride;
    ff = SigmaMatrix::Zero(matrix_size, matrix_size);
    if(use_fg)
    {   fg = SigmaMatrix::Zero(matrix_size, matrix_size);
        gf = SigmaMatrix::Zero(matrix_size, matrix_size);
    }
    if(use_gg)
        gg = SigmaMatrix::Zero(matrix_size, matrix_size);

    const double* R = lattice->R();
    const double* dR = lattice->dR();

    Rgrid.resize(matrix_size);
    dRgrid.resize(matrix_size);
    for(int i = 0; i < matrix_size; i++)
    {
        Rgrid[i] = R[start + i * stride];
        dRgrid[i] = dR[start + i * stride] * stride;
    }
}

unsigned int SigmaPotential::size() const
{   return (start + matrix_size * stride);
}

void SigmaPotential::IncludeLower(bool include_fg, bool include_gg)
{
    use_fg = include_fg;

    unsigned int new_matrix_size = use_fg? matrix_size: 0;
    fg = SigmaMatrix::Zero(new_matrix_size, new_matrix_size);
    gf = SigmaMatrix::Zero(new_matrix_size, new_matrix_size);

    use_gg = include_gg;

    new_matrix_size = use_gg? matrix_size: 0;
    gg = SigmaMatrix::Zero(new_matrix_size, new_matrix_size);
}

void SigmaPotential::AddToSigma(const SpinorFunction& s1, const SpinorFunction& s2, double coeff)
{
    // PRE: s1.size() & s2.size() >= size()
    // Map used subset of f1, f2 on to Eigen vectors
    EigenVectorMapped f1(s1.f.data()+start, matrix_size);
    EigenVectorMapped f2(s2.f.data()+start, matrix_size);

    ff.noalias() += coeff * f1 * f2.transpose();

    if(use_fg || use_gg)
    {
        EigenVectorMapped g1(s1.g.data()+start, matrix_size);
        EigenVectorMapped g2(s2.g.data()+start, matrix_size);

        if(use_fg)
        {   fg.noalias() += coeff * f1 * g2.transpose();
            gf.noalias() += coeff * g1 * f2.transpose();
        }
        if(use_gg)
            gg.noalias() += coeff * g1 * g2.transpose();
    }
}

void SigmaPotential::AddToSigma(const std::vector<double>& f1, const std::vector<double>& f2, double coeff)
{
    // PRE: s1.size() & s2.size() >= size()
    // Map used subset of f1, f2 on to Eigen vectors
    EigenVectorMapped mapped_f1(f1.data()+start, matrix_size);
    EigenVectorMapped mapped_f2(f2.data()+start, matrix_size);

    ff.noalias() += coeff * mapped_f1 * mapped_f2.transpose();
}

SpinorFunction SigmaPotential::ApplyTo(const SpinorFunction& a) const
{
    // PRE: a.size() >= size()
    SpinorFunction ret(a.Kappa(), start + matrix_size*stride);

    // Map used subset of a.f and a.g on to Eigen arrays (coefficient-wise multiplication).
    EigenArrayMapped fa(a.f.data()+start, matrix_size);
    EigenArrayMapped dr(lattice->dR()+start, matrix_size);

    // coefficient-wise multiplication
    Eigen::VectorXd fadr = (fa * dr * double(stride)).matrix();
    Eigen::VectorXd sigma_a_f = ff * fadr;

    if(use_fg)
    {
        EigenArrayMapped ga(a.g.data()+start, matrix_size);
        Eigen::VectorXd gadr = (ga * dr * double(stride)).matrix();

        // Add fg part to upper
        sigma_a_f += fg * gadr;

        // Add gf part to lower
        Eigen::VectorXd sigma_a_g = gf * fadr;

        // Add gg part to lower
        if(use_gg)
            sigma_a_g += gg * gadr;

        if(stride == 1)
            std::copy(sigma_a_g.data(), sigma_a_g.data()+matrix_size, ret.g.begin()+start);
        else
        {   const double* R = lattice->R();
            const std::vector<double> ag(sigma_a_g.data(), sigma_a_g.data()+matrix_size);
            Interpolator interp(Rgrid, dRgrid);
            double deriv;
            for(int i = start; i < ret.size(); i++)
                interp.Interpolate(ag, R[i], ret.g[i], deriv, 6);
        }
    }

    if(stride == 1)
        std::copy(sigma_a_f.data(), sigma_a_f.data()+matrix_size, ret.f.begin()+start);
    else
    {   const double* R = lattice->R();
        const std::vector<double> af(sigma_a_f.data(), sigma_a_f.data()+matrix_size);
        Interpolator interp(Rgrid, dRgrid);
        double deriv;
        for(int i = start; i < ret.size(); i++)
            interp.Interpolate(af, R[i], ret.f[i], deriv, 6);
    }

    return ret;
}

bool SigmaPotential::Read(const std::string& filename)
{
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
    {   return false;
    }

    fread(&matrix_size, sizeof(unsigned int), 1, fp);
    fread(&start, sizeof(unsigned int), 1, fp);

    fread(&use_fg, sizeof(bool), 1, fp);
    fread(&use_gg, sizeof(bool), 1, fp);

    resize_and_clear(start + matrix_size * stride);

    fread(ff.data(), sizeof(double), matrix_size * matrix_size, fp);

    if(use_fg)
    {   fread(fg.data(), sizeof(double), matrix_size * matrix_size, fp);
        fread(gf.data(), sizeof(double), matrix_size * matrix_size, fp);
    }

    if(use_gg)
        fread(gg.data(), sizeof(double), matrix_size * matrix_size, fp);

    fclose(fp);
    return true;
}

void SigmaPotential::Write(const std::string& filename) const
{
    FILE* fp = fopen(filename.c_str(), "wb");

    fwrite(&matrix_size, sizeof(unsigned int), 1, fp);
    fwrite(&start, sizeof(unsigned int), 1, fp);

    fwrite(&use_fg, sizeof(bool), 1, fp);
    fwrite(&use_gg, sizeof(bool), 1, fp);

    // Write data
    fwrite(ff.data(), sizeof(double), matrix_size * matrix_size, fp);

    if(use_fg)
    {   fwrite(fg.data(), sizeof(double), matrix_size * matrix_size, fp);
        fwrite(gf.data(), sizeof(double), matrix_size * matrix_size, fp);
    }
    if(use_gg)
        fwrite(gg.data(), sizeof(double), matrix_size * matrix_size, fp);

    fclose(fp);
}
