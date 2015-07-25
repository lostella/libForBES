/* 
 * File:   Quadratic.h
 * Author: Chung
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
#include <iostream>

/**
 * A Quadratic function of the form <code>Q(x) = 0.5 * x'*Q*x + q'*x</code>, where
 * <code>Q</code> is a square symmetric matrix and <code>q</code> is a vector. 
 * 
 * \ingroup Functions
 */
class Quadratic : public Function {
public:
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
    Quadratic(Matrix&Q);
    
    /**
     * Create a quadratic-plus-linear function of the form <code>f(x) = 0.5*x'*Q*x + q'*x</code>,
     * where <code>Q</code> is a square matrix and <code>q</code> is a vector.
     * @param Q A square matrix
     * @param q A vector
     */
    Quadratic(Matrix& Q, Matrix& q); // both Q and q    
    
    /**
     * The copy constructor of this class.
     * @param orig Instance of Quadratic to be copied.
     */
    Quadratic(const Quadratic& orig);
    
    /**
     * Destructor.
     */
    virtual ~Quadratic();

    /**
     * Category of this function. This is a quadratic function, therefore it 
     * returns <code>CAT_QUADRATIC</code>.
     * 
     * @return Category index.
     */
    int category();

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
     * status code which is equal to <code>STATUS_OK</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     * 
     */
    virtual int call(Matrix& x, double& f);
    
    virtual int callConj(const Matrix& x, double& f_star);
    
    virtual int callConj(const Matrix& x, double& f_star, Matrix& grad);



private:
    Matrix *Q;      /**< Matrix Q */
    Matrix *q;      /**< Vector q*/
    Matrix *L;      /**< Cholesky factor L */
    bool is_Q_eye;  /**< TRUE if Q is the identity matrix */
    bool is_q_zero; /**< TRUE is q is the zero vector */

    virtual int computeGradient(Matrix& x, Matrix& grad);


};

#endif	/* QUADRATIC_H */

