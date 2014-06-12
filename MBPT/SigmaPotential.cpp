#include "Include.h"
#include "SigmaPotential.h"

typedef Eigen::Map<const Eigen::VectorXd> EigenVectorMapped;
typedef Eigen::Map<const Eigen::ArrayXd> EigenArrayMapped;

SigmaPotential::SigmaPotential():
    start(0), matrix_size(0)
{}

SigmaPotential::SigmaPotential(unsigned int end_point, unsigned int start_point):
    start(start_point), matrix_size(0)
{
    resize_and_clear(end_point);
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
    matrix_size = new_size - start;
    ff = SigmaMatrix::Zero(matrix_size, matrix_size);
    if(use_fg)
    {   fg = SigmaMatrix::Zero(matrix_size, matrix_size);
        gf = SigmaMatrix::Zero(matrix_size, matrix_size);
    }
    if(use_gg)
        gg = SigmaMatrix::Zero(matrix_size, matrix_size);
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

SpinorFunction SigmaPotential::ApplyTo(const SpinorFunction& a, pLattice lattice) const
{
    // PRE: a.size() >= size()
    SpinorFunction ret(a.Kappa(), start+matrix_size);

    // Map used subset of a.f and a.g on to Eigen arrays (coefficient-wise multiplication).
    EigenArrayMapped fa(a.f.data()+start, matrix_size);
    EigenArrayMapped dr(lattice->dR()+start, matrix_size);

    // coefficient-wise multiplication
    Eigen::VectorXd fadr = (fa * dr).matrix();
    Eigen::VectorXd sigma_a_f = ff * fadr;

    if(use_fg)
    {
        EigenArrayMapped ga(a.g.data()+start, matrix_size);
        Eigen::VectorXd gadr = (ga * dr).matrix();

        // Add fg part to upper
        sigma_a_f += fg * gadr;

        // Add gf part to lower
        Eigen::VectorXd sigma_a_g = gf * fadr;

        // Add gg part to lower
        if(use_gg)
            sigma_a_g += gg * gadr;

        // Copy to ret
        std::copy(sigma_a_g.data(), sigma_a_g.data()+matrix_size, ret.g.begin()+start);
    }

    std::copy(sigma_a_f.data(), sigma_a_f.data()+matrix_size, ret.f.begin()+start);

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

    ff = SigmaMatrix::Zero(matrix_size, matrix_size);
    fread(ff.data(), sizeof(double), matrix_size * matrix_size, fp);

    if(use_fg)
    {   fg = SigmaMatrix::Zero(matrix_size, matrix_size);
        gf = SigmaMatrix::Zero(matrix_size, matrix_size);

        fread(fg.data(), sizeof(double), matrix_size * matrix_size, fp);
        fread(fg.data(), sizeof(double), matrix_size * matrix_size, fp);
    }
    else
    {   fg = SigmaMatrix::Zero(0, 0);
        gg = SigmaMatrix::Zero(0, 0);
    }

    if(use_gg)
    {   gf = SigmaMatrix::Zero(matrix_size, matrix_size);
        fread(fg.data(), sizeof(double), matrix_size * matrix_size, fp);
    }
    else
        gg = SigmaMatrix::Zero(0, 0);

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
