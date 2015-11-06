/* 
 * File:   LDLFactorization.h
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

#ifndef LDLFACTORIZATION_H
#define	LDLFACTORIZATION_H

#include <cstring>

#include "ForBESUtils.h"
#include "Matrix.h"
#include "FactoredSolver.h"
#include "ldl.h"

/**
 * 
 * \brief LDL factorization and solver.
 */
class LDLFactorization : public FactoredSolver {
public:

    /**
     * Creates an LDL factorizer given a %Matrix object.
     * 
     * @param m_matrix matrix to be factorized
     * 
     * \exception std::invalid_argument if an empty matrix is provided, or the 
     * matrix is of type <code>MATRIX_LOWERTR</code>, or if the matrix is not
     * square.
     * 
     * \exception std::logic_error the factorization cannot be applied to matrices
     * of type <code>MATRIX_DIAGONAL</code>. However, solving diagonal systems is
     * a trivial case.
     */
    explicit LDLFactorization(Matrix& m_matrix);

    /**
     * Destructor.
     */
    virtual ~LDLFactorization();

    virtual int factorize(void);

    /**
     * Solves a linear system using this instance of LDLFactorization.
     * 
     * \note Method #solve does not make use of the reference to the original matrix,
     * so it is not a problem if that matrix goes out of scope, is altered or deleted.
     * 
     * @param solution the solution of the linear system as an instance of <code>Matrix</code>
     * @param rhs The right-hand side vector or matrix
     * 
     * @return Returns <code>0</code> if the solution of the system has succeeded.
     * 
     * \sa FactoredSolver::solve
     */
    virtual int solve(const Matrix& rhs, Matrix& solution) const;

    double* getLDL() const;

    int* getIpiv() const;

private:


    double* LDL; /**< LDL factorization computed by lapack */
    int* ipiv;   /**< Pivots for the LDL factorization computed by lapack */

    /**
     * A sparse LDL factorization
     */
    typedef struct sparse_ldl_factor_struct {
        double * Lx;    /**< Values of L */
        int * Li;       /**< i-pointers of L */
        int * Lp;       /**< p-pointers of L */
        double * D;     /**< Diagonal part of the LDL factorization*/
    } sparse_ldl_factor;

    /**
     * Pointer to a sparse LDL factorization.
     */
    sparse_ldl_factor * m_sparse_ldl_factor;    
    

};

#endif	/* LDLFACTORIZATION_H */

