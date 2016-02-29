/* 
 * File:   Norm.h
 * Author: Pantelis Sopasakis
 *
 * Created on October 30, 2015, 5:32 PM
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

#ifndef NORM_H
#define	NORM_H

#include "Function.h"
#include <math.h>

/**
 * \class Norm
 * \brief Abstract norm class
 * \version version 0.1
 * \ingroup Functions
 * \date Created on October 30, 2015, 5:32 PM
 * \author Pantelis Sopasakis
 */
class Norm : public Function {
public:
    
    using Function::callConj;

    /**
     * The dual norm is defined as 
     * 
     * \f[
     *  \|x\|_\star = \sup_{\|z\|\leq 1} \langle z, x\rangle.
     * \f]
     * 
     * Knowing the dual norm allows us to define the conjugate function of a norm, that is
     * for \f$f(x)=\|x\|\f$, its conjugate is given by
     * 
     * \f[
     *  f^*(y) = 
     *  \begin{cases}
     *      0,&\text{if } \|y\|_\star \leq 1,\\
     *      \infty,&\text{otherwise}
     *  \end{cases}
     * \f]
     * 
     * The conjugate function of \f$f(x)=\|x\|\f$ is the indicator function of the
     * dual ball \f$\mathcal{B}_\star = \{y\mid \|y\|_\star \leq 1 \}\f$.
     * 
     * @param x point in \f$\mathbb{R}^n\f$ where the dual norm should be computed
     * @param norm the value of the dual norm at x, that is \f$\|x\|_\star\f$
     * @return 
     * status code which is equal to <code>STATUS_OK</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     */
    virtual int dualNorm(Matrix &x, double &norm);

    virtual int callConj(Matrix& x, double& f_star);

    /**
     * Norm category
     * @return 
     */
    virtual FunctionOntologicalClass category();

    /**
     * Default destructor of Norm
     */
    virtual ~Norm();
protected:
    Norm();

private:


};

#endif	/* NORM_H */

