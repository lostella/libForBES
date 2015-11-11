/* 
 * File:   CholeskyFactorization.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on August 4, 2015, 8:14 PM
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

#include "CholeskyFactorization.h"
#include "ForBESUtils.h"

CholeskyFactorization::CholeskyFactorization(Matrix& matrix) :
FactoredSolver(matrix) {
    m_L = NULL;
    m_factor = NULL;
    if (matrix.getNrows() != matrix.getNcols()){
        throw std::invalid_argument("CholeskyFactorization factorization can only be applied to square matrices");
    }
    if (matrix.getType() != Matrix::MATRIX_SPARSE) {
        this->m_L = new double[matrix.length()]();
    }
}

CholeskyFactorization::~CholeskyFactorization() {
    if (m_L != NULL) {
        delete[] m_L;
        m_L = NULL;
    }
    if (m_factor != NULL) {
        cholmod_free_factor(&m_factor, Matrix::cholmod_handle());
        m_factor = NULL;
    }
}

int CholeskyFactorization::factorize() {
    if (m_matrix_type == Matrix::MATRIX_SPARSE) {
        /* Cholesky decomposition of a SPARSE matrix: */
        if (m_matrix->m_sparse == NULL) {
            m_matrix->_createSparse();
        }
        /* analyze */
        m_factor = cholmod_analyze(m_matrix->m_sparse, Matrix::cholmod_handle());
        /* factorize */
        cholmod_factorize(m_matrix->m_sparse, m_factor, Matrix::cholmod_handle());
        /* Success: status = 0, else 1*/
        return (m_factor->minor == m_matrix->m_nrows) ? ForBESUtils::STATUS_OK : ForBESUtils::STATUS_NUMERICAL_PROBLEMS;
    } else { /* If this is any non-sparse matrix: */
        memcpy(m_L, m_matrix->getData(), m_matrix->length() * sizeof (double)); /* m_L := m_matrix.m_data */
        int info = ForBESUtils::STATUS_OK;
        if (m_matrix_type == Matrix::MATRIX_DENSE) { /* This is a dense matrix */
            info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', m_matrix_nrows, m_L, m_matrix_nrows);
#ifdef SET_L_OFFDIAG_TO_ZERO
            for (size_t i = 0; i < m_matrix_nrows; i++) {
                for (size_t j = i + 1; j < m_matrix_nrows; j++) {
                    L.set(i, j, 0.0);
                }
            }
#endif
        } else if (m_matrix_type == Matrix::MATRIX_SYMMETRIC) { /* This is a symmetric matrix */
            info = LAPACKE_dpptrf(LAPACK_COL_MAJOR, 'L', m_matrix_nrows, m_L);
        }
        return info;
    }
}

int CholeskyFactorization::solve(Matrix& rhs, Matrix& solution) const {
    if (m_matrix_type == Matrix::MATRIX_SPARSE) {
        cholmod_dense *x;

        /* cast the RHS as a cholmod_dense b = rhs */
        if (rhs.m_type == Matrix::MATRIX_DENSE) {

            cholmod_dense *b;
            b = cholmod_allocate_dense(rhs.m_nrows, rhs.m_ncols, rhs.m_nrows, CHOLMOD_REAL, Matrix::cholmod_handle());
            b->x = rhs.m_data;

            /* Solve - rhs is dense*/
            x = cholmod_solve(CHOLMOD_A, m_factor, b, Matrix::cholmod_handle());
            solution = Matrix(rhs.m_nrows, rhs.m_ncols);
            solution.m_delete_data = false;
            memcpy(solution.m_data, static_cast<double*> (x->x), rhs.m_nrows * rhs.m_ncols * sizeof (double));
            cholmod_free_dense(&x, Matrix::cholmod_handle());

        } else if (rhs.m_type == Matrix::MATRIX_SPARSE) {
            // still untested!
            cholmod_sparse * rhs_sparse;
            if (rhs.m_sparse == NULL) {
                const_cast<Matrix&> (rhs)._createSparse();
            }
            rhs_sparse = rhs.m_sparse;
            cholmod_sparse * result = cholmod_spsolve(CHOLMOD_LDLt, m_factor, rhs_sparse, Matrix::cholmod_handle());
            solution = Matrix(rhs.m_nrows, rhs.m_ncols, Matrix::MATRIX_SPARSE);
            solution.m_sparse = result;
            solution._createTriplet();
        } else {
            throw std::logic_error("Not supported");
        }
        return ForBESUtils::STATUS_OK;
    } else { /* the matrix to be factorized is not sparse */
        int info = ForBESUtils::STATUS_UNDEFINED_FUNCTION;
        solution = Matrix(rhs);
        if (m_matrix_type == Matrix::MATRIX_DENSE) {
            info = LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', m_matrix_nrows, rhs.m_ncols, m_L, m_matrix_nrows, solution.m_data, m_matrix_nrows);
        } else if (m_matrix_type == Matrix::MATRIX_SYMMETRIC) {
            info = LAPACKE_dpptrs(LAPACK_COL_MAJOR, 'L', m_matrix_nrows, rhs.m_ncols, m_L, solution.m_data, m_matrix_nrows);
        } else {
            throw std::invalid_argument("This matrix type is not supported - only DENSE, SPARSE and SYMMETRIC are supported");
        }
        return info;
    }
}


