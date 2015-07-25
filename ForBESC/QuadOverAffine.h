/* 
 * File:   QuadOverAffine.h
 * Author: Chung
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

/**
 * Quadratic-over-affine function.
 * 
 * This is the function <code>F(x) = x'*Q*x + q'*x + delta(x|Z)</code> where 
 * <code>delta(.|Z)</code> is the indicator function of an affine space <code>Z</code>
 * defined by <code>Z = {z: Az + b = 0}</code>.
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
    
    QuadOverAffine(const QuadOverAffine& orig);
    
    virtual ~QuadOverAffine();
    
    
private:
    
    QuadOverAffine();
    
    Matrix *Q;
    Matrix *q;
    Matrix *A;
    Matrix *b;
    
    Matrix *F;
};

#endif	/* QUADOVERAFFINE_H */

