/*
 * File:   IndBox.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 26, 2015, 5:22 PM
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

#ifndef INDBOX_H
#define	INDBOX_H

#include "Matrix.h"
#include "Function.h"
#include <algorithm>

/**
 * \class IndBox
 * \brief Indicator of a box
 * \version version 0.1
 * \ingroup Functions
 * \date Created on July 26, 2015, 5:22 PM
 * \author Pantelis Sopasakis
 * 
 * This class extends the class Function and implements the indicator function of
 * a box. 
 * 
 * A box is a set of the form \f$B_{[l,u]}=\{x \in \mathbb{R}^n: l\leq x \leq u\}\f$
 * and its indicator is defined as 
 * \f[
 * \delta (x\mid B_{[l,u]}) = \begin{cases}
 * 0, \text{ if } x\in B_{[l,u]},\\
 * \infty, \text{ otherwise}
 * \end{cases}
 * \f]
 * 
 * The proximal operator of this function is defined as 
 * \f[
 * \mathrm{prox}_{\gamma\delta(\cdot|B_{[l,u]})}(v) = \mathrm{mid}(v; l, u) = \min (\max(v, l), u)
 * \f]
 * These two functions are implemented and are available through #call and 
 * #callProx.
 * 
 * The conjugate of the indicator function of a box is the support function of
 * that box, that is
 * 
 * \f[
 * \begin{align}
 * \delta^*(y\mid B_{[l,u]}) &= \sup_{x\in B_{[l,u]}} \langle y, x\rangle = \sup_{l \leq x \leq u} \langle y, x\rangle\\
 * &= \sup_{l \leq x \leq u} \sum_{i=1}^{n}y_ix_i\\
 * &= \sum_{i=1}^{n} \sup_{l_i \leq x_i \leq u_i} y_ix_i\\
 * &= \sum_{i=1}^{n} \max \{ y_i l, y_i u \}
 * \end{align}
 * \f]
 * 
 * This function is available though #callConj.
 */
class IndBox : public Function {
public:

    using Function::call;
    using Function::callConj;

    /**
     * Constructor of instances of IndBox given a uniform lower bound and a uniform
     * upper bound. This then will be the indicator function of the set
     * \f[
     * B_{[l,u]} = \{x\in\mathbb{R}: l \leq x_i \leq u, \forall i=1,\ldots, n\},
     * \f]
     * where \f$l,u\in\mathbb{R}\f$ with \f$l\leq u\f$ are the uniform bounds
     * provided in this constructor.
     * @param uniform_lb Uniform lower bound
     * @param uniform_ub Uniform upper bound
     */
    IndBox(double& uniform_lb, double& uniform_ub);

    /**
     * Constructor of instances of IndBox given the lower and upper bounds as
     * instances of Matrix (column vectors). This then will be the indicator 
     * function of the set
     * \f[
     * B_{[l,u]} = \{x\in\mathbb{R}: l_i \leq x_i \leq u_i, \forall i=1,\ldots, n\},
     * \f]
     * where \f$l,u\in\mathbb{R}^n\f$ with \f$l_i\leq u_i\f$ are the bounds (column vectors,
     * instances of Matrix) provided in this constructor.
     * 
     * @param lb lower bound
     * @param ub upper bound
     */
    IndBox(Matrix& lb, Matrix& ub);

    /**
     * Default destructor
     */
    virtual ~IndBox();

    virtual int call(Matrix& x, double& f);

    virtual int callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    virtual int callProx(Matrix& x, double gamma, Matrix& prox);

    /**
     * Computes the conjugate of IndBox at a point <code>x</code> which is 
     * given by
     * 
     * \f[
     * \delta^*(y\mid B_{[l,u]}) = \sum_{i=1}^{n} \max \{ y_i l, y_i u \}
     * \f]
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
     */
    virtual int callConj(Matrix& x, double& f_star);

    virtual FunctionOntologicalClass category();


private:

    Matrix * m_lb;
    Matrix * m_ub;
    double * m_uniform_lb;
    double * m_uniform_ub;


};

#endif	/* INDBOX_H */

