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
 * \version 0.4
 * \author Pantelis Sopasakis
 * \ingroup LinSysSolver-group
 * \date November 11, 2015, 12:41 AM
 * 
 * A conjugate gradient solver which can be used to solver linear operator
 * equations of the form \f$T(x)=b\f$ for given \f$b\f$.
 * 
 * The conjugate gradient (CG) algorithm with a preconditioner \f$P\f$, which is
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
 * 7. Return \f$x\f$
 * 
 * Providing a preconditioner is optional. If no preconditioner is provided, it is 
 * assumed that \f$P\f$ is the identity operator, \f$P(x)=x\f$.
 * 
 * 
 * Systems of the form \f$Ax=b\f$, i.e., where \f$T(x)=Ax\f$ where \f$A\f$ is a 
 * Matrix can be solved using the linear operator MatrixOperator which wraps 
 * matrices as instances of LinearOperator.
 *  
 * Here is an example of use:
 * 
 * \code{.cpp}
 * size_t n = 20;
 * Matrix b = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 5.0);
 * Matrix A = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 0.1, Matrix::MATRIX_SYMMETRIC);
 * Matrix Y = MatrixFactory::MakeIdentity(n, 1.0);
 * A += Y;
 *
 * Matrix ID(n, n, Matrix::MATRIX_DIAGONAL);
 * for (size_t j = 0; j < n; ++j) {
 *     ID.set(j, j, 1 / A.get(j, j));
 * }
 * MatrixOperator Aop(A);
 * MatrixOperator M(ID);
 * 
 * size_t max_iter = n;
 * CGSolver solver(Aop, M, 1e-4, max_iter);
 * Matrix sol(n, 1);
 * int status = solver.solve(b, sol);
 * \endcode
 * 
 */
class CGSolver : public LinOpSolver {
public:

    /**
     * Constructs a new instance of CGSolver    
     * @param linop the underlying linear operator
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
     * @param linop linear operator which defined the system \f$T(x) = b\f$
     * @param preconditioner preconditioner as a linear operator
     * @param tolerance tolerance (default value, when other constructors are
     * used, is \f$10^{-4}\f$).
     * @param max_iterations maximum number of iterations after which the algorithm terminates
     * (the default value, when other constructors are used, is <code>500</code>).
     */
    CGSolver(LinearOperator& linop, LinearOperator& preconditioner, double tolerance, size_t max_iterations);


    /**
     * 
     * Solves the operator equation \f$T(x) = b\f$ for a given right-hand side \f$b\f$.
     * 
     * @param rhs the right-hand side of the equation
     * @param solution the solution to be computed
     * @return status code
     */
    virtual int solve(Matrix& rhs, Matrix& solution);

    int solve(Matrix& rhs, Matrix& solution, double tolerance, Matrix guess);

    /**
     * Default destructor.
     */
    virtual ~CGSolver();
    
    /**
     * Returns the infinity-norm of the last error recorded while the algorithm
     * was running. To assess the quality of the solution, it important to take into
     * account the exit flag of #solve as well as the number of iterations returned
     * by #last_num_iter.
     * @return last error
     * 
     * \sa #last_num_iter
     */
    double last_error() const;
    
    /**
     * The number of iterations of the algorithm on its last run.
     * @return number of iterations on last invocation.
     * 
     * \sa #last_num_iter
     */
    size_t last_num_iter() const;

private:

    LinearOperator * m_precond;
    double m_tolerance;
    double m_err;
    size_t m_max_iterations;
    size_t m_iterations_count;
    
    
    void init();   


protected:



};

#endif	/* CGSOLVER_H */

