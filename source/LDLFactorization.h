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

    explicit LDLFactorization(Matrix& m_matrix);

    virtual ~LDLFactorization();

    virtual int factorize(void);

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

