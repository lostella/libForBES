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

class Quadratic : public Function {
public:
    Quadratic();
    Quadratic(Matrix&); // provide just Q
    Quadratic(Matrix&,  Matrix&); // both Q and q    
    Quadratic(const Quadratic& orig);
    virtual ~Quadratic();

    int category();

    virtual int call( Matrix& x, float& f);
    virtual int callConj(const Matrix& x, float& f_star);
    virtual int callProx(const Matrix& x, float gamma, Matrix& prox, float f_at_prox);
    virtual int callProx(const Matrix& x, float gamma, Matrix& prox);   


protected:
    Matrix *Q; /**< Matrix Q */
    Matrix *q; /**< Vector q*/
    Matrix *L; /**< Cholesky factor L */
    bool is_Q_eye;  /**< TRUE if Q is the identity matrix */
    bool is_q_zero; /**< TRUE is q is the zero vector */
    
    virtual int computeGradient( Matrix& x, Matrix& grad);


};

#endif	/* QUADRATIC_H */

