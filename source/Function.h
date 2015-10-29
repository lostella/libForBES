/* 
 * File:   Function.h
 * Author: Pantelis Sopasakis
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
#include "ForBESUtils.h"
#include "FunctionOntologicalClass.h"
#include "FunctionOntologyRegistry.h"

/**
 * \class Function
 * \brief A ForBES function.
 * \version version 0.1
 * \ingroup Functions
 * \date Created on July 9, 2015, 3:35 AM
 * \author Pantelis Sopasakis
 * 
 * \todo Implement indBox, indPos, indNeg and indBall.
 * 
 * \todo Return, when possible, the Hessian as a linear operator.
 * 
 * This is a generic API for ForBES functions. Constructors for this class are
 * protected. Instances of this class can be created using any of its subclasses.
 */
class Function {
public:

    /**
     * Destructor for this class.
     */
    virtual ~Function();

    /**
     * Ontological categorization of the function.
     * 
     * @return function ontological class
     */
    virtual FunctionOntologicalClass category();

    /**
     * Returns the value of function f.
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
    virtual int call(Matrix& x, double& f);;

    /**
     * Same as <code>call(const Matrix& x, double& f)</code>, but this function returns
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
    virtual int call(Matrix& x, double& f, Matrix& grad); // returns also the gradient

    /**
     * Computes the proximal of this function at a point <code>x</code> with 
     * parameter <code>gamma</code>.
     * 
     * @param x The vector x where f(x) should be computed.
     * @param gamma The parameter <code>gamma</code> of <code>prox_(gamma*f)(x)</code>
     * @param prox The result of this operation
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     */
    virtual int callProx(const Matrix& x, double gamma, Matrix& prox);; // returns the value of prox_{gamma f}

    /**     
     * 
     * @param x The vector x where f(x) should be computed.
     * @param gamma The parameter <code>gamma</code> of <code>prox_(gamma*f)(x)</code>
     * @param prox The result of this operation
     * @param f_at_prox Value of this function at the proximal operator
     * 
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     */
    virtual int callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox);; // prox_{gamma f} and value-at-prox

    /**
     * Computes the conjugate of this function at a point <code>x</code>.
     * 
     * @param x The vector x where <code>f_star(x)</code> should be computed.
     * @param f_star the computed value <code>f_star(x)</code> 
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     */
    virtual int callConj(const Matrix& x, double& f_star);; // conjugate of f at x: f*(x)

    /**
     * Computes the conjugate of this function at a point <code>x</code> as well 
     * as the corresponding gradient.
     * 
     * @param x The vector x where <code>f_star(x)</code> should be computed.
     * @param f_star the computed value <code>f_star(x)</code> 
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     */
    virtual int callConj(const Matrix& x, double& f_star, Matrix& grad);; // Nabla f*(x)



private:



protected:
    /*
     * Constructors are protected.
     * It is not allowed to instantiate objects of this class directly.
     * Use its subclasses/implementations.
     */

    Function(); /**< Default constructor */

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
    virtual int computeGradient(Matrix& x, Matrix& grad);;

};

#endif	/* FUNCTION_H */

