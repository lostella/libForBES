/* 
 * File:   LDLFactorization_AAt.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on November 5, 2015, 1:12 AM
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

#include "S_LDLFactorization.h"
#include "ForBESUtils.h"
#include "LDLFactorization.h"

Matrix S_LDLFactorization::multiply_AAtr_betaI(Matrix& A, double beta) {
    size_t n = A.getNrows();
    size_t m = A.getNcols();
    Matrix result(n, n, Matrix::MATRIX_SYMMETRIC);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            result.set(i, j, 0.0);
            for (size_t k = 0; k < m; k++) {
                result.set(i, j, result.get(i, j) + A.get(i, k) * A.get(j, k));
            }
            if (i == j) {
                result.set(i, j, result.get(i, j) + beta);
            }
        }
    }
    return result;
}

S_LDLFactorization::S_LDLFactorization(Matrix& matrix, double beta) : FactoredSolver(matrix), m_beta(beta) {
    m_factor = NULL;
    m_delegated_solver = NULL;
}

int S_LDLFactorization::factorize() {
    if (m_matrix_type == Matrix::MATRIX_SPARSE) {
        double beta_temp[2];
        beta_temp[0] = m_beta;
        beta_temp[1] = 0.0;

        if (m_matrix.m_sparse == NULL) {
            m_matrix._createSparse();
        }

        m_matrix.m_sparse->stype = 0;
        m_factor = cholmod_analyze(m_matrix.m_sparse, Matrix::cholmod_handle());
        cholmod_factorize_p(m_matrix.m_sparse, beta_temp, NULL, 0, m_factor, Matrix::cholmod_handle());
        return (m_factor->minor == m_matrix.m_nrows) ? ForBESUtils::STATUS_OK : ForBESUtils::STATUS_NUMERICAL_PROBLEMS;
    } else if (m_matrix_type == Matrix::MATRIX_DENSE) {
        /* 
         * We here need to factorize a dense matrix 
         * We distinguish between two cases        
         * 1. A is short (more columns than rows)
         * 2. A is tall  (more rows than columns)
         */
        if (m_matrix.getNrows() <= m_matrix.getNcols()) {
            /* this is a ###SHORT### matrix */
            /*
             * Note: here we're creating matrix F on the fly and we're passing
             * its REFERENCE to LDLFactorization. Eventually, m_delegated_solver will
             * lose track of the original matrix.
             */
            Matrix F = multiply_AAtr_betaI(m_matrix, m_beta);
            m_delegated_solver = new LDLFactorization(F);
            int status = m_delegated_solver->factorize();
            return status;
        } else {
            /* this is a ~~~TALL~~~ matrix */
            m_matrix.transpose();
            /*
             * Note: here we're creating matrix F_tilde on the fly and we're passing
             * its REFERENCE to LDLFactorization. Eventually, m_delegated_solver will
             * lose track of the original matrix.
             */
            Matrix F_tilde = multiply_AAtr_betaI(m_matrix, m_beta);
            m_matrix.transpose();
            m_delegated_solver = new LDLFactorization(F_tilde);
            int status = m_delegated_solver->factorize();
            return status;
        }
    } else {
        throw std::invalid_argument("[uoe] Unsupported operation");
    }
}

int S_LDLFactorization::solve(const Matrix& rhs, Matrix& solution) const {
    solution = Matrix(rhs.m_nrows, rhs.m_ncols);
    if (m_matrix_type == Matrix::MATRIX_SPARSE) {
        if (m_factor == NULL) {
            throw std::invalid_argument(__FCT_MISS_EXCPT);
        }
        cholmod_dense *b;
        cholmod_dense *x;
        b = cholmod_allocate_dense(rhs.m_nrows, rhs.m_ncols, rhs.m_nrows, CHOLMOD_REAL, Matrix::cholmod_handle());
        b->x = rhs.m_data;
        x = cholmod_solve(CHOLMOD_A, m_factor, b, Matrix::cholmod_handle());
        solution.m_delete_data = false;
        memcpy(solution.m_data, static_cast<double*> (x->x), rhs.m_nrows * rhs.m_ncols * sizeof (double));
        cholmod_free_dense(&x, Matrix::cholmod_handle());
        return ForBESUtils::STATUS_OK;
    } else if (m_matrix_type == Matrix::MATRIX_DENSE) {
        if (m_delegated_solver == NULL) {
            throw std::invalid_argument(__FCT_MISS_EXCPT);
        }
        if (m_matrix_nrows <= m_matrix_ncols) {
            /* m_matrix is ###SHORT### and dense */
            int status = m_delegated_solver->solve(rhs, solution);
            return status;
        } else {
            /* m_matrix is ~~~TALL~~~ and dense */
            m_matrix.transpose();
            Matrix temp = m_matrix * const_cast<Matrix&>(rhs);
            m_matrix.transpose();
            Matrix c;
            int status = m_delegated_solver->solve(temp, c);
            if (status != ForBESUtils::STATUS_OK){
                return status;
            }
            solution = const_cast<Matrix&>(rhs) - m_matrix * c;
            double beta_inv = 1.0/m_beta;
            solution *= beta_inv;
            return ForBESUtils::STATUS_OK;
        }

    } else {
        throw std::invalid_argument("[uoe] Unsupported operation");
    }
}

S_LDLFactorization::~S_LDLFactorization() {
    if (m_factor != NULL) {
        cholmod_free_factor(&m_factor, Matrix::cholmod_handle());
        m_factor = NULL;
    }
    if (m_delegated_solver != NULL) {
        delete m_delegated_solver;
    }
}

