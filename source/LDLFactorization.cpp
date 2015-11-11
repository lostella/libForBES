/* 
 * File:   LDLFactorization.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on July 30, 2015, 3:02 AM
 * 
 * ForBES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *  
 * ForBES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with ForBES. If not, see <http://www.gnu.org/licenses/>.
 */

#include "LDLFactorization.h"

LDLFactorization::LDLFactorization(Matrix& matr) : FactoredSolver(matr) {
    this->LDL = NULL;
    this->ipiv = NULL;
    this->m_sparse_ldl_factor = NULL;
    this->m_matrix_type = m_matrix->getType();
    this->m_matrix_nrows = m_matrix->getNrows();
    if (matr.isEmpty()){
        throw std::invalid_argument("LDL factorization cannot be applied to empty matrices");
    }
    if (m_matrix_type == Matrix::MATRIX_LOWERTR) {
        throw std::invalid_argument("LDL factorization cannot be applied to non-symmetric matrices (e.g., Lower/Upper triangular)");
    }
    if (m_matrix_type != Matrix::MATRIX_DENSE &&
            matr.getType() != Matrix::MATRIX_SYMMETRIC &&
            matr.getType() != Matrix::MATRIX_SPARSE) { // Only DENSE and SYMMETRIC and SPARSE are supported!
        throw std::logic_error("This matrix type is not supported by LDLFactorization");
    }
    if (matr.getNrows() != matr.getNcols()) {
        throw std::invalid_argument("Matrix not square");
    }
    if (m_matrix_type == Matrix::MATRIX_SPARSE) {
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

int LDLFactorization::factorize() {
    int status = ForBESUtils::STATUS_UNDEFINED_FUNCTION;
    if (this->m_matrix_type == Matrix::MATRIX_DENSE) {
        status = LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'L', m_matrix_nrows, LDL, m_matrix_nrows, ipiv);
    } else if (this->m_matrix_type == Matrix::MATRIX_SYMMETRIC) {
        status = LAPACKE_dsptrf(LAPACK_COL_MAJOR, 'L', m_matrix_nrows, LDL, ipiv);
    } else if (this->m_matrix_type == Matrix::MATRIX_SPARSE) {
        // Factorize sparse matrix
        m_matrix->_createSparse();
        int * Parent = new int[m_matrix_nrows];
        int * Lnz = new int[m_matrix_nrows];
        int * Flag = new int[m_matrix_nrows];
        m_sparse_ldl_factor->Lp = new int[m_matrix_nrows + 1];
        ldl_symbolic(m_matrix_nrows,
                static_cast<int*>(m_matrix->m_sparse->p),
                static_cast<int*>(m_matrix->m_sparse->i),
                m_sparse_ldl_factor->Lp,
                Parent,
                Lnz,
                Flag, 
                NULL, NULL);
        int lnz = m_sparse_ldl_factor->Lp[m_matrix_nrows];
        int d;
        m_sparse_ldl_factor->Li = new int[lnz];
        m_sparse_ldl_factor->Lx = new double[lnz];
        m_sparse_ldl_factor->D = new double[m_matrix_nrows];

        double *Y = new double[m_matrix_nrows];
        int *Pattern = new int[m_matrix_nrows];

        d = ldl_numeric(m_matrix_nrows,
                static_cast<int*>(m_matrix->m_sparse->p),
                static_cast<int*>(m_matrix->m_sparse->i),
                static_cast<double*>(m_matrix->m_sparse->x),
                m_sparse_ldl_factor->Lp,
                Parent,
                Lnz,
                m_sparse_ldl_factor->Li,
                m_sparse_ldl_factor->Lx,
                m_sparse_ldl_factor->D,
                Y, 
                Pattern, 
                Flag, 
                NULL, NULL);
        
        delete[] Parent;
        delete[] Pattern;
        delete[] Y;
        delete[] Flag;
        delete[] Lnz;
        
        return d == m_matrix_nrows ? ForBESUtils::STATUS_OK : ForBESUtils::STATUS_NUMERICAL_PROBLEMS;

    } else {
        throw std::invalid_argument("This matrix type is not supported by LDLFactorization");
    }    
    return status;
}

int LDLFactorization::solve(Matrix& rhs, Matrix& solution) const {
    solution = Matrix(m_matrix_nrows, 1, Matrix::MATRIX_DENSE); // solution = rhs (DENSE)
    for (size_t i = 0; i < m_matrix_nrows; i++) {
        solution.set(i, 0, rhs.get(i, 0));
    }
    int status = ForBESUtils::STATUS_OK;
    if (Matrix::MATRIX_DENSE == this->m_matrix_type) {        
        status = LAPACKE_dsytrs(LAPACK_COL_MAJOR, 'L', m_matrix_nrows, 1, LDL, m_matrix_nrows, ipiv, solution.getData(), m_matrix_nrows);
    } else if (Matrix::MATRIX_SYMMETRIC == this->m_matrix_type) {
        status = LAPACKE_dsptrs(LAPACK_COL_MAJOR, 'L', m_matrix_nrows, 1, LDL, ipiv, solution.getData(), m_matrix_nrows);
    } else if (Matrix::MATRIX_SPARSE == this->m_matrix_type) {
        double * b = solution.getData();
        ldl_lsolve(m_matrix_nrows, b,
                m_sparse_ldl_factor->Lp,
                m_sparse_ldl_factor->Li,
                m_sparse_ldl_factor->Lx);
        ldl_dsolve(m_matrix_nrows, b, m_sparse_ldl_factor->D);
        ldl_ltsolve(m_matrix_nrows, b, m_sparse_ldl_factor->Lp, m_sparse_ldl_factor->Li, m_sparse_ldl_factor->Lx);
        status = ForBESUtils::STATUS_OK;
    }
    return status;
}

double* LDLFactorization::getLDL() const {
    return LDL;
}

int* LDLFactorization::getIpiv() const {
    return ipiv;
}