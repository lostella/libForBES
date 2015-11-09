/* 
 * File:   ConjugateFunction.h
 * Author: Pantelis Sopasakis
 *
 * Created on November 7, 2015, 3:15 AM
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

#ifndef CONJUGATEFUNCTION_H
#define	CONJUGATEFUNCTION_H

#include "Function.h"

/**
 * \class ConjugateFunction
 * \brief Conjugate of a given function
 * \version version 0.1
 * \ingroup Functions
 * \date Created on November 7, 2015, 3:15 AM
 * \author Pantelis Sopasakis
 * 
 * This class is a utility which provides access to the convex conjugate of a 
 * function which is passed to the constructor of ConjugateFunction by reference.
 * 
 * For a given function \f$s:X\to\mathbb{R}\cup \{+\infty\}\f$, this class 
 * provides access to \f$f(x) = s^*(x)\f$. Then,
 * 
 * 1. %Function \f$f(x) = s^*(x)\f$
 * 2. Gradient \f$\nabla f(x) = \nabla s^*(x)\f$
 * 3. Conjugate function \f$f^*(x) = s(x)\f$
 * 4. Conjugate gradient \f$\nabla f^*(x) = \nabla s(x)\f$
 * 5. Proximal operator \f$\mathrm{prox}_{\gamma s^*}(v) = v - \mathrm{prox}_{\gamma s}(v)\f$
 * 
 * See also the documentation of the various methods in this function
 * 
 */
class ConjugateFunction : public Function {
public:
    
    using Function::call;
    using Function::callConj;

    explicit ConjugateFunction(Function& funct);

    virtual ~ConjugateFunction();

    /**
     * For a given function \f$s\f$, provided in the constructor, this method 
     * computes its conjugate \f$s^*\f$.
     * @param x
     * @param f
     * @return 
     */
    virtual int call(Matrix& x, double& f);

    /**
     * For a given function \f$s\f$, provided in the constructor, this method 
     * computes its conjugate \f$s^*\f$ and its gradient \f$\nabla s^*\f$.
     * @param x
     * @param f
     * @return 
     */
    virtual int call(Matrix& x, double& f, Matrix& grad);

    /**
     * For a given function \f$s\f$, provided in the constructor, this method 
     * computes its value \f$s(x)=f^{**}(x)\f$ at a given point \f$x\in X\f$ as well
     * as the gradient \f$\nabla f^*(x) = \nabla s(x)\f$.
     * 
     * @param x The vector x where \f$f^*(x)\f$ should be computed.
     * @param f_star the computed value \f$f^*(x)\f$
     * @param grad the gradient of the conjugate function \f$\nabla f^*(x)\f$
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     * 
     * \exception std::invalid_argument an <code>invalid_argument</code> exception
     * is thrown in case the input %Matrix <code>x</code> or <code>grad</code> 
     * is of incompatible dimensions.
     */
    virtual int callConj(Matrix& x, double& f_star, Matrix& grad);

    /**
     * For a given function \f$s\f$, provided in the constructor, this method 
     * computes its value \f$s(x)\f$ at a given point \f$x\in X\f$.
     * 
     * @param x The vector x where \f$f^*(x)\f$ should be computed.
     * @param f_star the computed value \f$f^*(x)\f$
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     * 
     * \exception std::invalid_argument an <code>invalid_argument</code> exception
     * is thrown in case the input %Matrix <code>x</code> is of incompatible dimensions
     */
    virtual int callConj(Matrix& x, double& f_star);

    /**
     * For a given function \f$s\f$, provided in the constructor, this method 
     * computes the proximal gradient of its conjugate which is given by
     * \f[
     *  \mathrm{prox}_{\gamma s^*}(v) = v - \gamma \mathrm{prox}_{\gamma^{-1} s}(\gamma^{-1}v).
     * \f]
     *
     * @param x The vector x where \f$\mathrm{prox}_{\gamma f}(x)\f$ should be computed.
     * @param gamma The parameter \f$\gamma\f$ of \f$\mathrm{prox}_{\gamma f}\f$
     * @param prox The result of this operation
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     * 
     * \exception std::invalid_argument an <code>invalid_argument</code> exception
     * is thrown in case the function argument <code>x</code> and/or <code>prox</code>
     * are of incompatible dimensions.
     */
    virtual int callProx(Matrix& x, double gamma, Matrix& prox);

    /**
     * For a given function \f$s\f$, provided in the constructor, this method 
     * computes the proximal gradient of its conjugate which is given by
     * \f[
     *  \mathrm{prox}_{\gamma s^*}(v) = v - \gamma \mathrm{prox}_{\gamma^{-1} s}(\gamma^{-1}v).
     * \f]
     * and it also computes the value of \f$s^*\f$ at the proximal point, that is
     * \f$s^*(\mathrm{prox}_{\gamma s^*}(v))\f$.
     * 
     * @param x The vector x where \f$\mathrm{prox}_{\gamma f}(x)\f$ should be computed.
     * @param gamma The parameter \f$\gamma\f$ of \f$\mathrm{prox}_{\gamma f}\f$
     * @param prox The result of this operation
     * @param f_at_prox Value of this function at the proximal operator \f$f(\mathrm{prox}_{\gamma f}(x))\f$
     * 
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     * 
     * \exception std::invalid_argument an <code>invalid_argument</code> exception
     * is thrown in case the function argument <code>x</code> and/or <code>prox</code>
     * are of incompatible dimensions.
     */
    virtual int callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    /**
     * Ontological categorization of the function which is inferred from the 
     * category of the linked function \f$f\f$.
     * 
     * @return function ontological class
     */
    virtual FunctionOntologicalClass category();



private:

    Function& m_function;

};

#endif	/* CONJUGATEFUNCTION_H */

