/* 
 * File:   CholeskyFactorization.h
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

#ifndef CHOLESKYFACTORIZATION_H
#define	CHOLESKYFACTORIZATION_H

#include "FactoredSolver.h"
#include <cstring>

/**
 * 
 * \brief Cholesky factorization and solver.
 */
class CholeskyFactorization : public FactoredSolver {
public:

    CholeskyFactorization(Matrix& m_matrix) :
    FactoredSolver(m_matrix) {
        if (m_matrix.getType() != Matrix::MATRIX_SPARSE) {
            this->m_L = new double[m_matrix.length()]();
        }
    }

    virtual ~CholeskyFactorization();

    /**
     * Computes the Cholesky factorization of this matrix. 
     * 
     * <p>If applied on a <code>MATRIX_DENSE</code>
     * matrix, then it is assumed it is symmetric (but there is no verification) and
     * only its lower triangular part is considered. Notice that the Cholesky factorization
     * can only be applied to symmetric and positive definite matrices.</p>  
     *      
     * @return status code. Returns <code>0</code> if the factorization succeeded.
     *      
     */
    virtual int factorize();

    /**
     * Solves the linear system <code>L*L'*x = b</code>, where <code>L</code> is the
     * current factorization (implicitly), <code>b</code> is a given right-hand side vector or matrix
     * (given as an instance of <code>Martrix</code>) and <code>x</code> is the
     * solution which is returned by this method.
     * 
     * Note that if this is a <code>MATRIX_DENSE</code> matrix, then it is assumed it is 
     * lower triangular (but there is no verification) and only its lower triangular 
     * part is considered.
     * 
     * @param solution the solution of the linear system as an instance of <code>Matrix</code>
     * @param rhs The right-hand side vector or matrix
     * 
     * @return Returns <code>0</code> if the solution of the system has succeeded.
     */
    virtual int solve(const Matrix& rhs, Matrix& solution) const;



private:
    double * m_L = NULL;
    cholmod_factor * m_factor = NULL;

};

#endif	/* CHOLESKYFACTORIZATION_H */

