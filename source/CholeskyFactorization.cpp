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

CholeskyFactorization::CholeskyFactorization(Matrix& m_matrix) :
FactoredSolver(m_matrix) {
    m_L = NULL;
    if (m_matrix.getType() != Matrix::MATRIX_SPARSE) {
        this->m_L = new double[m_matrix.length()]();
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
    size_t n = m_matrix.getNrows();
    if (m_matrix.getType() == Matrix::MATRIX_SPARSE) {
        /* Cholesky decomposition of a SPARSE matrix: */
        if (m_matrix.m_sparse == NULL) {
            m_matrix._createSparse();
        }
        m_factor = cholmod_analyze(m_matrix.m_sparse, Matrix::cholmod_handle()); // analyze
        cholmod_factorize(m_matrix.m_sparse, m_factor, Matrix::cholmod_handle()); // factorize       
        return (m_factor->minor == m_matrix.m_nrows) ? 0 : 1; /* Success: status = 0, else 1*/
    } else { /* If this is any non-sparse matrix: */
        memcpy(m_L, m_matrix.getData(), m_matrix.length() * sizeof (double)); /* m_L := m_matrix.m_data */
        int info = ForBESUtils::STATUS_OK;
        if (m_matrix.getType() == Matrix::MATRIX_DENSE) { /* This is a dense matrix */
            info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', n, m_L, n);
#ifdef SET_L_OFFDIAG_TO_ZERO
            for (size_t i = 0; i < n; i++) {
                for (size_t j = i + 1; j < n; j++) {
                    L.set(i, j, 0.0);
                }
            }
#endif
        } else if (m_matrix.getType() == Matrix::MATRIX_SYMMETRIC) { /* This is a symmetric matrix */
            info = LAPACKE_dpptrf(LAPACK_COL_MAJOR, 'L', n, m_L);
        }
        return info;
    }
}

int CholeskyFactorization::solve(const Matrix& rhs, Matrix& solution) const {
    if (m_matrix.getType() == Matrix::MATRIX_SPARSE) {
        cholmod_dense *x;
        cholmod_dense *b;
        b = cholmod_allocate_dense(rhs.m_nrows, rhs.m_ncols, rhs.m_nrows, CHOLMOD_REAL, Matrix::cholmod_handle());
        for (size_t k = 0; k < rhs.getNrows(); k++) {
            ((double*) b->x)[k] = rhs.get(k, 0);
        }
        x = cholmod_solve(CHOLMOD_A, m_factor, b, Matrix::cholmod_handle());
        solution = Matrix(rhs.m_nrows, rhs.m_ncols);
        for (size_t k = 0; k < x->nzmax; k++) {
            solution.m_data[k] = ((double*) x->x)[k];
        }
        cholmod_free_dense(&x, Matrix::cholmod_handle());
        return 0;
    } else {
        int info = ForBESUtils::STATUS_UNDEFINED_FUNCTION;
        solution = Matrix(rhs);
        size_t n = m_matrix.getNrows();
        if (m_matrix.getType() == Matrix::MATRIX_DENSE) {
            info = LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', n, rhs.m_ncols, m_L, n, solution.m_data, n);
        } else if (m_matrix.getType() == Matrix::MATRIX_SYMMETRIC) {
            info = LAPACKE_dpptrs(LAPACK_COL_MAJOR, 'L', n, rhs.m_ncols, m_L, solution.m_data, n);
        }
        return info;
    }
}


