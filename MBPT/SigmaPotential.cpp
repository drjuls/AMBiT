#include "Include.h"
#include "SigmaPotential.h"
#include "Universal/Interpolator.h"

namespace Ambit
{
typedef Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<>> EigenVectorMapped;
typedef Eigen::Map<const Eigen::ArrayXd, Eigen::Unaligned, Eigen::InnerStride<>> EigenArrayMapped;

SigmaPotential::SigmaPotential(pLattice lattice):
    LatticeObserver(lattice), start(0), matrix_size(0), stride(4)
{}

SigmaPotential::SigmaPotential(pLattice lattice, unsigned int end_point, unsigned int start_point, unsigned int stride):
    LatticeObserver(lattice), start(start_point), matrix_size(0), stride(stride)
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

void SigmaPotential::Alert()
{
    if(lattice->size() < size())
    {
        matrix_size = (lattice->size() - start)/stride;
        ff.noalias() = ff.topLeftCorner(matrix_size, matrix_size);
        if(use_fg)
        {   fg.noalias() = fg.topLeftCorner(matrix_size, matrix_size);
            gf.noalias() = gf.topLeftCorner(matrix_size, matrix_size);
        }
        if(use_gg)
            gg.noalias() = gg.topLeftCorner(matrix_size, matrix_size);

        Rgrid.resize(matrix_size);
        dRgrid.resize(matrix_size);
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
    EigenVectorMapped f1(s1.f.data()+start, matrix_size, Eigen::InnerStride<>(stride));
    EigenVectorMapped f2(s2.f.data()+start, matrix_size, Eigen::InnerStride<>(stride));

    ff.noalias() += coeff * f1 * f2.transpose();

    if(use_fg || use_gg)
    {
        EigenVectorMapped g1(s1.g.data()+start, matrix_size, Eigen::InnerStride<>(stride));
        EigenVectorMapped g2(s2.g.data()+start, matrix_size, Eigen::InnerStride<>(stride));

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
    EigenVectorMapped mapped_f1(f1.data()+start, matrix_size, Eigen::InnerStride<>(stride));
    EigenVectorMapped mapped_f2(f2.data()+start, matrix_size, Eigen::InnerStride<>(stride));

    ff.noalias() += coeff * mapped_f1 * mapped_f2.transpose();
}

SpinorFunction SigmaPotential::ApplyTo(const SpinorFunction& a) const
{
    // PRE: a.size() >= size()
    SpinorFunction ret(a.Kappa(), start + matrix_size*stride);

    // Map used subset of a.f and a.g on to Eigen arrays (coefficient-wise multiplication).
    EigenArrayMapped fa(a.f.data()+start, matrix_size, Eigen::InnerStride<>(stride));
    EigenArrayMapped dr(lattice->dR()+start, matrix_size, Eigen::InnerStride<>(stride));

    // coefficient-wise multiplication
    Eigen::VectorXd fadr = (fa * dr * double(stride)).matrix();
    Eigen::VectorXd sigma_a_f = ff * fadr;

    if(use_fg)
    {
        EigenArrayMapped ga(a.g.data()+start, matrix_size, Eigen::InnerStride<>(stride));
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
    FILE* fp = file_err_handler->fopen(filename.c_str(), "rb");
    if(!fp)
    {   return false;
    }

    file_err_handler->fread(&matrix_size, sizeof(unsigned int), 1, fp);
    file_err_handler->fread(&start, sizeof(unsigned int), 1, fp);
    file_err_handler->fread(&stride, sizeof(unsigned int), 1, fp);

    file_err_handler->fread(&use_fg, sizeof(bool), 1, fp);
    file_err_handler->fread(&use_gg, sizeof(bool), 1, fp);

    resize_and_clear(start + matrix_size * stride);

    file_err_handler->fread(ff.data(), sizeof(double), matrix_size * matrix_size, fp);

    if(use_fg)
    {   file_err_handler->fread(fg.data(), sizeof(double), matrix_size * matrix_size, fp);
        file_err_handler->fread(gf.data(), sizeof(double), matrix_size * matrix_size, fp);
    }

    if(use_gg)
        file_err_handler->fread(gg.data(), sizeof(double), matrix_size * matrix_size, fp);

    file_err_handler->fclose(fp);
    return true;
}

void SigmaPotential::Write(const std::string& filename) const
{
    if(ProcessorRank == 0)
    {
        FILE* fp = file_err_handler->fopen(filename.c_str(), "wb");

        file_err_handler->fwrite(&matrix_size, sizeof(unsigned int), 1, fp);
        file_err_handler->fwrite(&start, sizeof(unsigned int), 1, fp);
        file_err_handler->fwrite(&stride, sizeof(unsigned int), 1, fp);

        file_err_handler->fwrite(&use_fg, sizeof(bool), 1, fp);
        file_err_handler->fwrite(&use_gg, sizeof(bool), 1, fp);

        // Write data
        file_err_handler->fwrite(ff.data(), sizeof(double), matrix_size * matrix_size, fp);

        if(use_fg)
        {   file_err_handler->fwrite(fg.data(), sizeof(double), matrix_size * matrix_size, fp);
            file_err_handler->fwrite(gf.data(), sizeof(double), matrix_size * matrix_size, fp);
        }
        if(use_gg)
            file_err_handler->fwrite(gg.data(), sizeof(double), matrix_size * matrix_size, fp);

        file_err_handler->fclose(fp);
    }
}
}
