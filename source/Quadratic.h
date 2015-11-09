/* 
 * File:   Quadratic.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 9, 2015, 3:36 AM
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

#ifndef QUADRATIC_H
#define	QUADRATIC_H

#include "Function.h"
#include "Matrix.h"
#include "CholeskyFactorization.h"
#include <iostream>

/**
 * 
 * \class Quadratic
 * \brief A quadratic function <code>F(x) = 0.5 * x'*Q*x + q'*x</code>
 * \version 0.1
 * \author Pantelis Sopasakis
 * \date Created on July 24, 2015, 4:55 PM
 * A Quadratic function of the form <code>Q(x) = 0.5 * x'*Q*x + q'*x</code>, where
 * <code>Q</code> is a square symmetric matrix and <code>q</code> is a vector. 
 * 
 * \ingroup Functions
 * 
 * A quadratic function is a function of the form
 * 
 * \f[
 * F(x) = \frac{1}{2}x'Qx + q'x,
 * \f]
 * 
 * where <code>Q</code> is a symmetric positive definite matrix whose conjugate is given by
 * 
 * \f[
 * F^*(x^*) = (x^*-q)'Q^{-1}(x^*-q).
 * \f]
 * 
 * Here is a simple example:
 * \code
 * size_t n = 10;
 * Matrix Q = MatrixFactory::MakeRandomSparse(n, n, 20, 0.0, 1.0);
 * Matrix Eye = MatrixFactory::MakeIdentity(n, 10.0);
 * 
 * Matrix Qt(Q);        // Qt= Q
 * Qt.transpose();      // Qt = Qt'
 * Q += Qt;             // Q = Q + Qt (this will make Q a symmetric matrix)
 * 
 * Q += Eye;            // with this we ensure Q is positive definite
 *
 * Function *F = new Quadratic(Q); // F(x) = 0.5*x'Qx
 * Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
 * 
 * double f;            // value F(x)
 * double f_star;       // value F*(x) (conjugate of F at x)
 * Matrix grad;         // gradient of F* and x
 * 
 * int status = F->call(x, f);
 * status = F->callConj(x, f_star, grad);
 * 
 * std::cout << grad;   // print out the gradient
 * \endcode
 * 
 * The invocation of <code>callConj</code> involves the computation of a Cholesky
 * factor of <code>Q</code> which is stored internally in the instance of our 
 * quadratic function.
 */
class Quadratic : public Function {
public:
    
    using Function::call;
    using Function::callConj;
    
    /**
     * Create a trivial quadratic function with zero Hessian and
     * zero linear term.
     */
    Quadratic();

    /**
     * Create a quadratic function of the form <code>f(x) = 0.5*x'*Q*x</code>, where
     * <code>Q</code> is a square matrix.
     * 
     * @param Q Matrix Q.
     */
    explicit Quadratic(Matrix&Q);

    /**
     * Create a quadratic-plus-linear function of the form <code>f(x) = 0.5*x'*Q*x + q'*x</code>,
     * where <code>Q</code> is a square matrix and <code>q</code> is a vector.
     * @param Q A square matrix
     * @param q A vector
     */
    Quadratic(Matrix& Q, Matrix& q); // both Q and q    


    /**
     * Destructor.
     */
    virtual ~Quadratic();
    

    virtual FunctionOntologicalClass category();

  

    void setQ(Matrix& Q);
    
    void setq(Matrix& q);

    /**
     * Returns the value of function f which is computed as <code>Q(x)=0.5*x'*Q*x + q'*x</code>.
     * 
     * @param x The vector x where f(x) should be computed.
     * 
     * @param f The computed value of f(x)
     * 
     * @return
     * Status code which is equal to <code>STATUS_OK</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     * 
     */
    virtual int call(Matrix& x, double& f);

    virtual int call(Matrix& x, double& f, Matrix& grad);
    
    virtual int call(Matrix& x, double& f, Matrix& grad, Matrix& hessian);

    virtual int callConj(Matrix& x, double& f_star);

    virtual int callConj(Matrix& x, double& f_star, Matrix& grad);    


private:
    Matrix *m_Q;                  /**< Matrix Q */
    Matrix *m_q;                  /**< Vector q*/
    FactoredSolver *m_solver;     /**< Cholesky factor L */
    bool m_is_Q_eye;              /**< TRUE if Q is the identity matrix */
    bool m_is_q_zero;             /**< TRUE is q is the zero vector */

    /**
     * Computes the gradient of this function at a given vector x. 
     * @param x The vector x where the gradient of f should be computed.
     * @param grad the computed gradient at x
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     */
    virtual int computeGradient(Matrix& x, Matrix& grad);


};

#endif	/* QUADRATIC_H */

