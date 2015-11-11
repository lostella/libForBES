/* 
 * File:   FactoredSolver.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 30, 2015, 3:18 AM
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

#ifndef FACTOREDSOLVER_H
#define	FACTOREDSOLVER_H

#include "Matrix.h"
#include "LinSysSolver.h"
#include "MatrixSolver.h"

/**
 * \class FactoredSolver
 * \brief An abstract factorization-based solver for linear systems.
 * \version version 0.1
 * \date Created on July 30, 2015, 3:18 AM
 * \author Pantelis Sopasakis
 * \ingroup LinSysSolver-group
 * 
 * FactoredSolver is a solver for linear systems of the form \f$Ax=b\f$ where a
 * factorization of matrix \f$A\f$ is used. This class exports its functionality
 * using two methods: #factorize and #solve. Objects of this class are instantiated
 * provided the matrix \f$A\f$ for which a reference is stored inside the object.
 * 
 * \sa LinSysSolver
 */
class FactoredSolver : public MatrixSolver {
public:

    /**
     * Creates a new instance of FactoredSolver given a reference to the matrix
     * to be factorized.
     * @param m_matrix new instance of FactoredSolver.
     */
    explicit FactoredSolver(Matrix& matrix);

    virtual ~FactoredSolver();

    /**
     * An abstract method which factorizes matrix \f$A\f$.
     * @return status code. The method returns \link ForBESUtils::STATUS_OK STATUS_OK\endlink 
     * if the invocation has succeeded or 
     * \link ForBESUtils::STATUS_NUMERICAL_PROBLEMS STATUS_NUMERICAL_PROBLEMS\endlink
     * if numerical problems have hindered the factorization of the matrix.
     * 
     * \sa #solve
     */
    virtual int factorize(void) = 0;

    /**
     * Solves the linear system using the system's matrix factorization 
     * and returns its solution and a status code.
     * 
     * \pre Always call #factorize before you call solve.
     * 
     * 
     * @param rhs the right-hand side of the linear equation
     * @param solution the solution of the linear system which is computed using the 
     * matrix factorization.
     * @return status code. The method returns \link ForBESUtils::STATUS_OK STATUS_OK\endlink 
     * if the invocation has succeeded or 
     * \link ForBESUtils::STATUS_NUMERICAL_PROBLEMS STATUS_NUMERICAL_PROBLEMS\endlink
     * if numerical problems have hindered the computation of a solution.
     * 
     * \sa #factorize
     */
    virtual int solve(Matrix& rhs, Matrix& solution) const = 0;

private:


};

#endif	/* FACTOREDSOLVER_H */

