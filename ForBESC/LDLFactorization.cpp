/* 
 * File:   LDLFactorization.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on July 30, 2015, 3:02 AM
 */

#include "LDLFactorization.h"


LDLFactorization::LDLFactorization(Matrix& matr) : FactoredSolver(matr) {
    if (matr.getType() == Matrix::MATRIX_LOWERTR) {
        throw std::invalid_argument("LDL factorization cannot be applied to non-symmetric matrices (e.g., Lower/Upper triangular)");
    }
    if (matr.getType() != Matrix::MATRIX_DENSE &&
            matr.getType() != Matrix::MATRIX_SYMMETRIC &&
            matr.getType() != Matrix::MATRIX_SPARSE) { // Only DENSE and SYMMETRIC and SPARSE are supported!
        throw std::logic_error("This matrix type is not supported by LDLFactorization");
    }
    if (matr.getNrows() != matr.getNcols()) {
        throw std::invalid_argument("Matrix not square");
    }
    if (matr.getType() == Matrix::MATRIX_SPARSE) {
        m_sparse_ldl_factor = new sparse_ldl_factor;
        return;
    }
    this->LDL = new double[matr.length()];
    this->ipiv = new int[matr.getNrows()];
    memcpy(this->LDL, matr.getData(), matr.length() * sizeof (double));
}

LDLFactorization::~LDLFactorization() {
    if (this->LDL != NULL) {
        delete[] LDL;
    }
    if (this->ipiv != NULL) {
        delete[] this->ipiv;
    }
}

int LDLFactorization::solve(const Matrix& rhs, Matrix& solution) const {
    size_t n = m_matrix.getNrows();
    solution = Matrix(rhs);
    int status;
    if (Matrix::MATRIX_DENSE == m_matrix.getType()) {
        status = LAPACKE_dsytrs(LAPACK_COL_MAJOR, 'L', n, 1, LDL, n, ipiv, solution.getData(), n);
    } else if (Matrix::MATRIX_SYMMETRIC == m_matrix.getType()) {
        status = LAPACKE_dsptrs(LAPACK_COL_MAJOR, 'L', n, 1, LDL, ipiv, solution.getData(), n);
    }
    return status;
}

double* LDLFactorization::getLDL() const {
    return LDL;
}

int* LDLFactorization::getIpiv() const {
    return ipiv;
}

int LDLFactorization::factorize() {
    size_t n = m_matrix.getNrows();
    int status = ForBESUtils::STATUS_UNDEFINED_FUNCTION;
    if (m_matrix.getType() == Matrix::MATRIX_DENSE) {
        status = LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'L', n, LDL, n, ipiv);
    } else if (m_matrix.getType() == Matrix::MATRIX_SYMMETRIC) {
        status = LAPACKE_dsptrf(LAPACK_COL_MAJOR, 'L', n, LDL, ipiv);
    } else if (m_matrix.getType() == Matrix::MATRIX_SPARSE) {
        // Factorize sparse matrix
        m_matrix._createSparse();
        int n = m_matrix.getNrows();
        int Parent[n];
        int Lnz[n];
        int Flag[n];
        m_sparse_ldl_factor->Lp = new int[n + 1];
        ldl_symbolic(n, 
                (int*)(m_matrix.m_sparse->p), 
                (int*)(m_matrix.m_sparse->i), 
                m_sparse_ldl_factor->Lp, 
                Parent, 
                Lnz, 
                Flag, NULL, NULL);
        int llz = m_sparse_ldl_factor->Lp[n];
        std::cout << "llz = " << llz << "\n";
    }
    return status;
}

