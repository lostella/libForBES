/* 
 * File:   Norm1.h
 * Author: Pantelis Sopasakis
 *
 * Created on October 30, 2015, 4:05 PM
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

#ifndef NORM1_H
#define	NORM1_H

#include "Function.h"
#include "Norm.h"

/**
 * \class Norm1
 * \brief %Norm-1 loss function
 * \version version 0.1
 * \ingroup Functions
 * \date Created on October 30, 2015, 4:05 PM
 * \author Pantelis Sopasakis
 * 
 * This class extends the class Function and implements the norm-1
 * function which is a function \f$g:\mathbb{R}^n\to\mathbb{R}_+\f$ defined as 
 * 
 * \f[
 *  g(x) = \mu\|x\|_1,
 * \f]
 *
 * where
 * 
 * \f[
 *  \|x\|_1 = \sum_{i=1}^{n} |x_i|.
 * \f]
 * 
 * The proximal operator of norm-1 is given by
 * 
 * \f[
 *  \mathrm{prox}_{\gamma \mu \|\cdot\|_1}(v)_i = 
 *  \begin{cases}
 *      v_i - \gamma \mu,   &\text{if } v_i \geq \gamma\mu,\\
 *      0,                  &\text{if } |v_i| \leq \gamma\mu,\\
 *      v_i + \gamma\mu,    &\text{if } v_i \leq - \gamma\mu
 *  \end{cases}
 * \f]
 */
class Norm1 : public Norm {
    
public:
    
    /**
     * Construct a new instance of Norm1 with \f$\mu=1\f$.
     */
    Norm1();
    
    /**
     * Construct a new instance of Norm1 with a given parameter \f$\mu>0\f$.
     * @param mu Parameter \f$\mu\f$
     */
    Norm1(double mu);
    
    /**
     * Destructor for Norm1.
     */
    virtual ~Norm1();
    

    /**
     * This norm is defined as \f$\|x\|=\mu\|x\|_1\f$, so its dual is by definition
     * \f[
     *  \|y\|_* = \sup_{\mu\|z\|_1\leq 1} \langle z, y\rangle = \frac{1}{\mu}\|x\|_\infty,
     * \f]
     * given that the dual of norm-1 is norm-infinity.
     * 
     * @param x point in \f$\mathbb{R}^n\f$ where the dual norm should be computed
     * @param norm the value of the dual norm at x, that is \f$\|x\|_*\f$
     * @return 
     * status code which is equal to <code>STATUS_OK</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     */
    virtual int dualNorm(const Matrix& x, double& norm);

    virtual int call(Matrix& x, double& f);    

    virtual int callProx(const Matrix& x, double gamma, Matrix& prox);
    
    virtual int callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox);
    
    virtual FunctionOntologicalClass category();


private:
    double m_mu;
};

#endif	/* NORM1_H */

