/* 
 * File:   Function.h
 * Author: chung
 *
 * Created on July 9, 2015, 3:35 AM
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

#ifndef FUNCTION_H
#define	FUNCTION_H

#include "Matrix.h"

class Function {
public:
    /**
     * Method has succeeded.
     */
    const static int STATUS_OK;
    /**
     * Method is undefined.
     */
    const static int STATUS_UNDEFINED_FUNCTION;
    /**
     * The result is unreliable, or could not be computed because
     * of numerical errors.
     */
    const static int STATUS_NUMERICAL_PROBLEMS;        
    
    /**
     * A quadratic function.
     */
    const static int CAT_QUADRATIC;
    
    
    
    
    Function();
    Function(const Function& orig);
    virtual ~Function();
    
    /**
     * The function category. Functions may define a function category so that
     * the caller can know what type of function this is. For instance quadratic
     * functions (e.g., our implementation <code>Quadratic</code>) return 
     * <code>Function::CAT_QUADRATIC</code>. Users can define and use their own
     * category indices. 
     * 
     * @return Function category as <code>int</code>
     */
    virtual int category() = 0; // abstract method
    
    /**
     * Returns the value of function f.
     * 
     * @param x The vector x where f(x) should be computed.
     * 
     * @param f The computed value of f(x)
     * 
     * @return
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     * 
     */
    virtual int call(const Matrix& x, float& f) const =0; 
    
    /**
     * Same as <code>call(const Matrix& x, float& f)</code>, but this function returns
     * also the gradient of f at x.
     * 
     * @param x The vector x where f(x) should be computed.
     * 
     * @param f The computed value of f(x)
     * 
     * @param grad The gradient of f at x.
     * 
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     */
    virtual int call(const Matrix& x, float& f, Matrix& grad) const; // returns also the gradient
    
    /**
     * 
     * @param x
     * @param gamma
     * @param prox
     * @return 
     */
    virtual int callProx(const Matrix& x, float gamma, Matrix& prox) =0; // returns the value of prox_{gamma f}
       
    virtual int callProx(const Matrix& x, float gamma, Matrix& prox, float f_at_prox) =0; // prox_{gamma f} and value-at-prox
    
    virtual int callConj(const Matrix& x, float& f_star) =0; // conjugate of f at x: f*(x)


private:

protected:
    virtual int computeGradient(const Matrix& x, Matrix& grad) const =0; 
    
    //virtual int computeFunctionAtProx(const Matrix& x, float gamma, const Matrix& prox, float f_at_prox) const =0;

};

#endif	/* FUNCTION_H */

