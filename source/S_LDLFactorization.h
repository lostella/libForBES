/* 
 * File:   LDLFactorization_AAt.h
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

#ifndef LDLFACTORIZATION_AAT_H
#define	LDLFACTORIZATION_AAT_H

#include "FactoredSolver.h"
#include "ForBESUtils.h"
#include "LDLFactorization.h"

#define __FCT_MISS_EXCPT "[factorization_missing] It seems #factorize has not been invoked, or it failed"

/**
 * 
 * \class S_LDLFactorization
 * \brief LDL factorization of A'A+bI
 * \version version 0.1
 * \date Created on November 5, 2015, 1:12 AM
 * \author Pantelis Sopasakis
 * \ingroup LinSysSolver-group
 * 
 * Let \f$A\in\mathbb{R}^{n\times m}\f$ be a non necessarily square matrix. This 
 * class, which is derived from the abstract class FactoredSolver, decomposes the
 * matrix \f$F=AA'+\beta I\f$ for a given scalar \f$\beta\f$ without the need to 
 * pre-compute \f$F\f$. 
 * 
 * In case \f$A\f$ is a short matrix (it has more columns than rows), then we compute
 * \f$AA'+\beta I\f$ and we delegate its factorization to LDLFactorization.
 * 
 * If \f$A\f$ is a tall matrix (it has more rows than columns) then, using the 
 * <a href="https://en.wikipedia.org/wiki/Woodbury_matrix_identity">Woodbury matrix identity</a> 
 * we have
 * 
 * \f[
 *  (AA'+\beta I)^{-1} = \beta^{-1}I - \beta^{-1}A(\beta I + A'A)^{-1}A'.
 * \f]
 * 
 * Notice that matrix \f$\tilde{F}=\beta I + A'A\f$ is of smaller size than \f$F\f$.
 * The above can be concisely written as 
 * 
 * \f[
 * F^{-1} = \beta^{-1}(I-A\tilde{F}^{-1}A')
 * \f]
 * 
 * Then, to compute a \f$z\in\mathbb{R}^{n}\f$ which solves
 * 
 * \f[
 *  Fz = t,
 * \f]
 * 
 * we simply do \f$z=\beta^{-1}(I-A\tilde{F}^{-1}A')t\f$ or
 * 
 * \f[
 *  z = \beta^{-1}(t+Ac),
 * \f]
 * 
 * where \f$c\f$ is the solution of the following system
 * \f[
 *  \tilde{F} c  = A't,
 * \f]
 * which is determined using the LDL factorization of \f$\tilde{F}\f$.
 * 
 * This class is powered by <a href="http://faculty.cse.tamu.edu/davis/suitesparse.html">SuiteSparse</a> 
 * for sparse matrices.
 */
class S_LDLFactorization : public FactoredSolver {
public:

    /**
     * Creates a new instance of S_LDLFactorization.
     * 
     * @param matrix Any matrix of type <code>MATRIX_DENSE</code> or <code>MATRIX_SPARSE</code>
     * @param beta A positive scalar
     */
    S_LDLFactorization(Matrix& matrix, double beta);

    virtual ~S_LDLFactorization();

    virtual int factorize();

    /**
     * Computes the solution of the linear system
     * 
     * \f[
     *  (AA' + \beta I) x = b,
     * \f]
     * 
     * for a given \f$b\f$.

     * 
     * \exception std::invalid_argument If the RHS is of incompatible dimensions, or
     * it is neither of \link Matrix::MATRIX_SPARSE sparse\endlink nor
     * \link Matrix::MATRIX_DENSE dense\endlink type.
     * 
     * @param rhs the right-hand side of the linear equation \f$(AA'+ \beta I) x = b\f$
     * 
     * @param solution the solution of the linear system \f$(AA'+\beta I) x = b\f$ 
     * which is computed using this matrix factorization.
     * 
     * @return status code. The method returns \link ForBESUtils::STATUS_OK STATUS_OK\endlink 
     * if the invocation has succeeded or 
     * \link ForBESUtils::STATUS_NUMERICAL_PROBLEMS STATUS_NUMERICAL_PROBLEMS\endlink
     * if numerical problems have hindered the computation of a solution.
     * 
     * \attention In case the matrix to be factorized is \link Matrix::MATRIX_DENSE dense\endlink
     * and <em>tall</em> (i.e., it has more rows than columns), it is imperative 
     * that the reference passed in the constructor of 
     * S_LDLFactorization remains available when #solve is called.
     * \par 
     * So, for instance, the following code would throw an error:
     * \code{.cpp}
     * // create a tall dense matrix:
     * size_t n = 500;
     * size_t m = 10;
     * Matrix * A = new Matrix(n, m, Matrix::MATRIX_DENSE);
     * 
     * // construct a S_LDLFactorization object
     * double beta = 0.690;
     * FactoredSolver * solver = new S_LDLFactorization(A, beta);
     * int status;
     * status = solver -> factorize();
     * 
     * // delete A
     * delete A; 
     * Matrix b(n,1); 
     * Matrix x;
     * 
     * // The solver cannot access A
     * status = solver -> solve(b, x);
     * \endcode
     * Note that %S_LDLFactorization does not keep a hard copy of the matrix
     * object that is passed to it, but just a reference.
     * \par
     * There are certain more subtle cases when this may happen such as the following
     * \code{.cpp}
     * class MyClass {
     *  public:
     *      MyClass(){
     *          // create a tall dense matrix:
     *          size_t n = 500;
     *          size_t m = 10;
     *          // Matrix A has the scope of this constructor
     *          Matrix A(n, m, Matrix::MATRIX_DENSE);
     *          double beta = 1.3523;
     *          // m_solver will hold a reference to A which in turn which 
     *          // will go out of scope after the constructor returns
     *          m_solver = new S_LDLFactorization(A, beta);
     *          int status = m_solver -> factorize();
     *      }
     * 
     *      void myMethod(){
     *          Matrix b = ...;
     *          Matrix x = ...;
     *          m_solver -> solve(b, x);
     *      }
     * 
     *  private:
     *      FactoredSolver * m_solver;
     * }
     * \endcode
     * A solution to the above pathological case is to store matrix objects in a
     * way that remain accessible from the %S_LDLFactorization object when #solve
     * is invoked.
     */
    virtual int solve(Matrix& rhs, Matrix& solution);

private:

    /**
     * A cholmod_factor used when m_matrix is sparse
     */
    cholmod_factor * m_factor;
    /**
     * Scalar beta
     */
    double m_beta;
    
    /**
     * A delegated LDL solver (used when m_matrix is dense)
     */
    FactoredSolver * m_delegated_solver;

    /**
     * Performs AA'+beta*I for dense matrices. The result will be a 
     * symmetric matrix (type <code>MATRIX_SYMMETRIC</code>).
     * 
     * @param A given matrix
     * @param beta scalar beta
     * @return matrix AA'
     */
    static Matrix multiply_AAtr_betaI(Matrix& A, double beta);

};

#endif	/* LDLFACTORIZATION_AAT_H */

