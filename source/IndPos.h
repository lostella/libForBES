/* 
 * File:   IndPos.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 26, 2015, 4:36 PM
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

#ifndef INDPOS_H
#define	INDPOS_H

#include "Function.h"
#include "IndBox.h"


/**
 * \class IndPos
 * \brief %Indicator of the positive orthant
 * \version 0.1
 * \author Pantelis Sopasakis
 * \date Created on July 26, 2015, 4:36 PM
 * 
 * This class implements the indicator function of a translated positive orthant,
 * that is the set
 *
 * \f[
 *    K_+^l = \{x\in\mathbb{R}^n, x\geq l\}
 * \f]
 *
 * Then, the indicator function of \f$K^l_+\f$ is defined as
 * 
 * \f[
 * \delta(x\mid K^l_+) = \begin{cases}
 * 0, & \text{ if } x\in K_+^l,\\
 * +\infty, & \text{ otherwise}
 * \end{cases}
 * \f]
 * 
 * The proximal operator of this indicator function is then
 * 
 * \f[
 * \mathrm{prox}_{\gamma \delta(\cdot\mid K^l_+)}(v) = \max(v,l),
 * \f]
 * 
 * and the conjugate function is
 * 
 * \f[
 * \begin{align}
 * \delta^*(y\mid K^l_+) &= \sup_{x\geq l} \langle x, y \rangle\\
 *                       &= \sup_{x\geq l} \sum_{i=1}^{n} y_i x_i\\
 *                       &= \sum_{i=1}^{n} \sup_{x_i\geq l_i} y_i x_i,
 * \end{align}
 * \f]
 * 
 * where
 * 
 * \f[
 * \sup_{x_i\geq l_i} y_i x_i = \begin{cases}
 * +\infty, &\text{ if } y_i > 0,\\
 * y_i l_i,   &\text{ otherwise}
 * \end{cases}
 * \f]
 * 
 * These functions are implemented as #call, #callProx and #callConj respectively.
 * 
 * 
 * \ingroup Functions
 */
class IndPos : public Function {
public:
    using Function::call;
    using Function::callConj;
    
    /**
     * Constructor for instances of IndPos given a uniform lower bound. This will
     * then be the indicator of the set
     * \f[
     * K_+^l = \{x\in\mathbb{R}^n: x_i \geq l, \forall i=1,\ldots, n\},
     * \f]
     * where \f$l\in\mathbb{R}\f$ is the given uniform lower bound.
     * @param uniform_lb uniform lower bound
     */
    explicit IndPos(double& uniform_lb);
    
    /**
     * Constructor for an instance of IndPos with a given lower bound as as instance
     * of Matrix (column vector). This will then be the indicator of the set
     * \f[
     * K_+^l = \{x\in\mathbb{R}^n: x_i \geq l_i, \forall i=1,\ldots, n\},
     * \f]
     * where \f$l\in\mathbb{R}^n\f$ is the given lower bound.
     * @param lb lower bound
     */
    explicit IndPos(Matrix& lb);
    
    /**
     * Constructor for a new instance of IndPos assuming that the lower 
     * bound is the zero vector, therefore, this becomes the indicator of the
     * positive cone:
     * 
     * \f[
     *  K_+ = \{x\in\mathbb{R}^n, x\geq 0\}
     * \f]
     * 
     */
    IndPos();
               
    /**
     * Default destructor
     */
    virtual ~IndPos();
    
    virtual int call(Matrix& x, double& f);

    virtual int callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    virtual int callProx(Matrix& x, double gamma, Matrix& prox);
    
    virtual int callConj(Matrix& y, double& f_star);
    
    virtual FunctionOntologicalClass category();



private:
    Matrix * m_lb;
    double * m_uniform_lb; 
    
    

};

#endif	/* INDPOS_H */

