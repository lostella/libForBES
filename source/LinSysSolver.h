/* 
 * File:   LinSysSolver.h
 * Author: Pantelis Sopasakis
 *
 * Created on November 7, 2015, 3:47 PM
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

#ifndef LINSYSSOLVER_H
#define	LINSYSSOLVER_H

#include "Matrix.h"
#include "LinearOperator.h"
#include "MatrixFactory.h"

/**
 * \class LinSysSolver
 * \brief An abstract solver for linear systems
 * \version version 0.0
 * \date Created on November 7, 2015, 3:47 PM
 * \author Pantelis Sopasakis
 * \ingroup LinSysSolver-group
 * 
 * LinSysSolver is a solver interface for linear systems of the form \f$Ax=b\f$,
 * where \f$A\f$ is a matrix, or systems of the form \f$T(x)=b\f$, where \f$T\f$
 * is a linear operator (instance of LinearOperator).
 */
class LinSysSolver {
public:
        
    
    virtual ~LinSysSolver();

    /**
     * Solves the linear system.
     * 
     * @param rhs the right-hand side of the linear equation
     * @param solution the solution of the linear system which is computed using the 
     * matrix factorization.
     * @return status code. The method returns \link ForBESUtils::STATUS_OK STATUS_OK\endlink 
     * if the invocation has succeeded or 
     * \link ForBESUtils::STATUS_NUMERICAL_PROBLEMS STATUS_NUMERICAL_PROBLEMS\endlink
     * if numerical problems have hindered the computation of a solution.
     * 
     */
    virtual int solve(Matrix& rhs, Matrix& solution) const = 0;

private:

protected:
    
    LinSysSolver();
    
       

};

#endif	/* LINSYSSOLVER_H */

