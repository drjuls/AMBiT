#ifndef HAMILTONIAN_MATRIX_H
#define HAMILTONIAN_MATRIX_H

#include "RelativisticConfiguration.h"
#include "NonRelConfiguration.h"
#include "HartreeFock/HFOperator.h"
#include "Level.h"
#include "Universal/Enums.h"
#include "Universal/Matrix.h"
#include "ManyBodyOperator.h"
#include "MBPT/OneElectronIntegrals.h"
#include "MBPT/TwoElectronCoulombOperator.h"
#include "MBPT/Sigma3Calculator.h"
#include <Eigen/Eigen>

namespace Ambit
{
typedef ManyBodyOperator<pHFIntegrals, pTwoElectronCoulombOperator> TwoBodyHamiltonianOperator;
typedef std::shared_ptr<TwoBodyHamiltonianOperator> pTwoBodyHamiltonianOperator;

typedef ManyBodyOperator<pHFIntegrals, pTwoElectronCoulombOperator, pSigma3Calculator> ThreeBodyHamiltonianOperator;
typedef std::shared_ptr<ThreeBodyHamiltonianOperator> pThreeBodyHamiltonianOperator;

/** The dimensions of HamiltonianMatrix is set by the RelativisticConfigList.
    It is generally size N * N, where N = relconfigs->NumCSFs(), however it also supports a "non-square" matrix
    with dimensions (Nsmall, N) where Nsmall = relconfigs->NumCSFsSmall(). In this case it stores a trapezoid,
    rather than a triangle, because the lower right (N-Nsmall) square is zero (except on the diagonal).
 */
class HamiltonianMatrix : public Matrix
{
public:
    /** Hamiltonian with one and two-body operators. */
    HamiltonianMatrix(pHFIntegrals hf, pTwoElectronCoulombOperator coulomb, pRelativisticConfigList relconfigs);
    /** Hamiltonian with additional effective three-body operator. */
    HamiltonianMatrix(pHFIntegrals hf, pTwoElectronCoulombOperator coulomb, pSigma3Calculator sigma3, pConfigListConst leadconfigs, pRelativisticConfigList relconfigs);
    virtual ~HamiltonianMatrix();

    virtual void MatrixMultiply(int m, double* b, double* c) const;
    virtual void GetDiagonal(double* diag) const;

    /** Generate Hamiltonian matrix. */
    virtual void GenerateMatrix(unsigned int configs_per_chunk = 4);

    /** Print upper triangular part of matrix (text). Lower triangular part is zeroed. */
    friend std::ostream& operator<<(std::ostream& stream, const HamiltonianMatrix& matrix);

    /** Read binary file if it exists. Return success. */
    virtual bool Read(const std::string& filename);

    /** Write binary file. */
    virtual void Write(const std::string& filename) const;

    /** Return proportion of elements that have magnitude greater than epsilon. */
    virtual double PollMatrix(double epsilon = 1.e-15) const;

    /** Solve the matrix that has been generated. Note that this may destroy the matrix. */
    virtual LevelVector SolveMatrix(pHamiltonianID hID, unsigned int num_solutions);

#ifdef AMBIT_USE_SCALAPACK
    /** Solve using ScaLAPACK. Note that this destroys the matrix.
        If use_energy_limit is true, only return eigenvalues below energy_limit (up to num_solutions of them).
     */
    virtual LevelVector SolveMatrixScalapack(pHamiltonianID hID, unsigned int num_solutions, bool use_energy_limit, double energy_limit = 0.0);

    virtual LevelVector SolveMatrixScalapack(pHamiltonianID hID, double energy_limit)
    {   return SolveMatrixScalapack(hID, N, true, energy_limit);
    }
#endif

    /** Clear matrix and recover memory. */
    virtual void Clear() { chunks.clear(); }

protected:
    pRelativisticConfigList configs;
    pTwoBodyHamiltonianOperator H_two_body;

    pConfigListConst leading_configs;           //!< Leading configs for sigma3
    pThreeBodyHamiltonianOperator H_three_body; //!< Three-body operator is null if sigma3 not used

    unsigned int Nsmall;            //!< For non-square CI, the smaller matrix size

protected:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;

    /** MatrixChunk is a rectangular section of the lower triangular part of the HamiltonianMatrix.
        The top left corner of the section is at (start_row, 0).
        The number of rows is num_rows, and the section goes to the diagonal of the Hamiltonian matrix (or Nsmall if chunk is in the extra part).
        The rows correspond to a number of RelativisticConfigurations in configs, as distributed by GenerateMatrix().
        RelativisticConfigurations included are [config_indices.first, config_indices.second).
        The matrix section "diagonal" is a square on the diagonal of the Hamiltonian that is outside Nsmall.
     */
    class MatrixChunk
    {
    public:
        MatrixChunk(unsigned int config_index_start, unsigned int config_index_end, unsigned int row_start, unsigned int num_rows, unsigned int Nsmall, bool big_chunk):
            start_row(row_start), num_rows(num_rows), is_big_chunk(big_chunk)
        {
            config_indices.first = config_index_start;
            config_indices.second = config_index_end;
            chunk = RowMajorMatrix::Zero(num_rows, mmin(start_row + num_rows, Nsmall));

            if(Nsmall < start_row + num_rows)
            {   unsigned int diagonal_size = mmin(num_rows, start_row + num_rows - Nsmall);
                diagonal = RowMajorMatrix::Zero(diagonal_size, diagonal_size);
            }
        }

        std::pair<unsigned int, unsigned int> config_indices;
        unsigned int start_row;
        unsigned int num_rows;
        RowMajorMatrix chunk;
        RowMajorMatrix diagonal;
        bool is_big_chunk;

        /** Make upper triangle part of the matrix chunk match the lower. */
        void Symmetrize()
        {
            if(start_row < chunk.cols())
            {
                for(unsigned int i = 0; i < chunk.cols() - start_row - 1; i++)
                    for(unsigned int j = i+start_row+1; j < chunk.cols(); j++)
                        chunk(i, j) = chunk(j - start_row, i + start_row);
            }

            if(diagonal.size())
            {
                for(unsigned int i = 0; i < diagonal.rows() - 1; i++)
                    for(unsigned int j = i + 1; j < diagonal.cols(); j++)
                        diagonal(i, j) = diagonal(j, i);
            }
        }
    };

    std::vector<MatrixChunk> chunks;
    unsigned int most_chunk_rows;
};

}
#endif
