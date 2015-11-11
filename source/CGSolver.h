/* 
 * File:   CGSolver.h
 * Author: Pantelis Sopasakis
 *
 * Created on November 11, 2015, 12:41 AM
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

#ifndef CGSOLVER_H
#define	CGSOLVER_H

#include "LinSysSolver.h"
#include "MatrixOperator.h"
#include "LinOpSolver.h"

/**
 * \class CGSolver
 * \brief Conjugate gradient solver
 * \version 0.0
 * \author Pantelis Sopasakis
 * \ingroup LinSysSolver-group
 * \date November 11, 2015, 12:41 AM
 * 
 * A conjugate gradient solver which can be used to solver linear operator
 * equations of the form \f$T(x)=b\f$ for given \f$b\f$.
 * 
 * The conjugate gradient (CG) algorithm with a predonditioner \f$P\f$, which is
 * itself a linear operator, is defined by the following iteration:
 * 
 * 1. \f$x\leftarrow x_0\f$
 * 2. \f$r\leftarrow b-T(x)\f$
 * 3. \f$z\leftarrow P(r)\f$
 * 4. \f$p\leftarrow z\f$
 * 5. Do:
 *      1. \f$\alpha \leftarrow \langle r, z\rangle/\langle p, T(p)\rangle\f$
 *      2. \f$x\leftarrow x + \alpha p\f$
 *      3. \f$r_{+} \leftarrow r - \alpha T(p)\f$
 *      4. If \f$\|r\|<\varepsilon\f$, break
 *      5. \f$z_{+} \leftarrow P(r_{+})\f$
 *      6. \f$\beta = \langle z_{+}, r_{+}\rangle/\langle z, r\rangle\f$
 *      7. \f$p \leftarrow z_{+} + \beta p\f$
 *      8. \f$z \leftarrow z_{+}\f$
 *      9. \f$r \leftarrow r_{+}\f$
 * 6. Repeat until the desired accuracy or the maximum number of iterations is reached.
 * 
 * 
 */
class CGSolver : public LinOpSolver {
public:

    /**
     * Constructs a new instance of CGSolver
     * @param linop
     */
    explicit CGSolver(LinearOperator& linop);

    /**
     * 
     * @param linop linear operator which defined the system \f$T(x) = b\f$
     * @param preconditioner preconditioner as a linear operator
     */
    CGSolver(LinearOperator& linop, LinearOperator& preconditioner);


    /**
     * 
     * Solves the operator equation \f$T(x) = b\f$ for a given right-hand side \f$b\f$.
     * 
     * @param rhs
     * @param solution
     * @return 
     */
    virtual int solve(Matrix& rhs, Matrix& solution) const;

    int solve(Matrix& rhs, Matrix& solution, double tolerance, Matrix guess);

    virtual ~CGSolver();

private:

    LinearOperator * m_precond;


protected:



};

#endif	/* CGSOLVER_H */

