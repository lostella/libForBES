/* 
 * File:   QuadOverAffine.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 24, 2015, 4:55 PM
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

#ifndef QUADOVERAFFINE_H
#define	QUADOVERAFFINE_H

#include "Function.h"
#include "FactoredSolver.h"

/**
 * \class QuadOverAffine
 * \brief %Function <code>F(x) = 0.5*x'*Q*x + q'*x + delta(x|Z)</code>
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on July 24, 2015, 4:55 PM
 * 
 * Quadratic-over-affine function.
 * 
 * This is the function \f$F(x) = \frac{1}{2}x'Qx + q'x + \delta(x|Z)\f$ where 
 * \f$\delta(\cdot|Z)\f$ is the indicator function of an affine space \f$Z\f$
 * defined by \f$Z = \{z: Az - b = 0\}\f$
 * 
 * This class implements only the computation of the conjugate of quadratic-over-affine
 * as well as its gradient.
 * 
 * \ingroup Functions
 */
class QuadOverAffine : public Function {
public:
    
    /**
     * Define a new quadratic-over-affine function.
     * @param Q
     * @param q
     * @param A
     * @param b
     */
    QuadOverAffine(Matrix& Q, Matrix& q, Matrix& A, Matrix& b);    
    
    virtual ~QuadOverAffine();
    
    virtual int callConj(const Matrix& y, double& f_star, Matrix& grad);

    
private:
    
    QuadOverAffine();
    
    Matrix *Q = NULL;   /**< Matrix Q (Hessian) */
    Matrix *q = NULL;   /**< Vector q (Linear term) */
    Matrix *A = NULL;   /**< Matrix A */
    Matrix *b = NULL;   /**< Matrix b */
    
    Matrix *F = NULL;                   /**< Matrix <code>F = [Q A'; A 0]</code> */
    Matrix *sigma = NULL;
    FactoredSolver * Fsolver = NULL;    /**< Factorizer for matrix F */
};

#endif	/* QUADOVERAFFINE_H */

