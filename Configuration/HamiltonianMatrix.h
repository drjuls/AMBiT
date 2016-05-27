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

typedef ManyBodyOperator<pHFIntegrals, pTwoElectronCoulombOperator> TwoBodyHamiltonianOperator;
typedef std::shared_ptr<TwoBodyHamiltonianOperator> pTwoBodyHamiltonianOperator;

typedef ManyBodyOperator<pHFIntegrals, pTwoElectronCoulombOperator, pSigma3Calculator> ThreeBodyHamiltonianOperator;
typedef std::shared_ptr<ThreeBodyHamiltonianOperator> pThreeBodyHamiltonianOperator;

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
    virtual void GenerateMatrix();

    /** Print upper triangular part of matrix (text). Lower triangular part is zeroed. */
    friend std::ostream& operator<<(std::ostream& stream, const HamiltonianMatrix& matrix);

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

protected:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;

    /** MatrixChunk is a rectangular section of the HamiltonianMatrix.
        The top left corner of the section is on the diagonal at (start_row, start_row).
        The number of rows is num_rows, and the section goes to column N (the right edge of the Hamiltonian matrix).
        The rows correspond to a number of RelativisticConfigurations in configs, as distributed by GenerateMatrix().
        RelativisticConfigurations included are [config_indices.first, config_indices.second).
     */
    class MatrixChunk
    {
    public:
        MatrixChunk(unsigned int config_index_start, unsigned int config_index_end, unsigned int row_start, unsigned int num_rows, unsigned int N):
            start_row(row_start), num_rows(num_rows)
        {   config_indices.first = config_index_start;
            config_indices.second = config_index_end;
            chunk = RowMajorMatrix::Zero(num_rows, N-start_row);
        }

        std::pair<unsigned int, unsigned int> config_indices;
        unsigned int start_row;
        unsigned int num_rows;
        RowMajorMatrix chunk;

        /** Make lower triangle part of the matrix chunk match the upper. */
        void Symmetrize()
        {
            for(unsigned int i = 1; i < num_rows; i++)
                for(unsigned int j = 0; j < i; j++)
                    chunk(i, j) = chunk(j, i);
        }
    };

    std::vector<MatrixChunk> chunks;
    unsigned int most_chunk_rows;
};

#endif
