/* 
 * File:   LDLFactorization_AAt.cpp
 * Author: chung
 * 
 * Created on November 5, 2015, 1:12 AM
 */

#include "S_LDLFactorization.h"
#include "ForBESUtils.h"

S_LDLFactorization::S_LDLFactorization(Matrix& matrix, double beta) : FactoredSolver(matrix), m_beta(beta) {
    m_factor = NULL;

}

int S_LDLFactorization::factorize() {
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
}

int S_LDLFactorization::solve(const Matrix& rhs, Matrix& solution) const {
    cholmod_dense *b;
    cholmod_dense *x;
    b = cholmod_allocate_dense(rhs.m_nrows, rhs.m_ncols, rhs.m_nrows, CHOLMOD_REAL, Matrix::cholmod_handle());
    b->x = rhs.m_data;
    x = cholmod_solve(CHOLMOD_A, m_factor, b, Matrix::cholmod_handle());
    solution = Matrix(rhs.m_nrows, rhs.m_ncols);
    solution.m_delete_data = false;
    memcpy(solution.m_data, static_cast<double*> (x->x), rhs.m_nrows * rhs.m_ncols * sizeof (double));
    cholmod_free_dense(&x, Matrix::cholmod_handle());
    return ForBESUtils::STATUS_OK;
}

S_LDLFactorization::~S_LDLFactorization() {
}

