/* 
 * File:   Norm2.h
 * Author: Pantelis Sopasakis
 *
 * Created on October 30, 2015, 6:54 PM
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

#ifndef NORM2_H
#define	NORM2_H

#include "Norm.h"
#include <math.h>

/**
 * \class Norm2
 * \brief Euclidean norm
 * \version version 0.1
 * \ingroup Functions
 * \date Created on October 30, 2015, 6:54 PM
 * \author Pantelis Sopasakis
 * 
 * This class defines the function
 * 
 * \f[
 * f(x) = \mu\|x\|_2,
 * \f]
 * 
 * for vectors \f$x\in\mathbb{R}^n\f$ and scalars \f$\mu>0\f$. %Function \f$f\f$ is
 * a norm \f$\|x\| = f(x)\f$ for all \f$\mu>0\f$ with dual
 * 
 * \f[
 *  \|x\|_* = \frac{1}{\mu}\|x\|_2.
 * \f]
 * 
 * This allows us to define the conjugate function of this function which is
 * the indicator of the \f$\|\cdot\|_*\f$-ball, that is
 * 
 *  \f[
 *  f^*(y) = 
 *  \begin{cases}
 *      0,&\text{if } \|y\|_\star \leq 1,\\
 *      \infty,&\text{otherwise}
 *  \end{cases}
 * \f]
 * 
 * The proximal operator of this function is given by
 * 
 * \f[
 * \mathrm{prox}_{\gamma f}(v) = \begin{cases}
 * 0, &\text{ if } \|v\| \leq \gamma \mu,\\
 * (1-\frac{\gamma\mu}{\|v\|})v,&\text{ otherwise}
 * \end{cases}
 * \f]
 */
class Norm2 : public Norm {
public:

    using Function::call;
    using Norm::callConj;

    /**
     * Default constructor. Here, it is assumed that \f$\mu=1\f$. It 
     * constructs an instance of the function \f$f(x) = \mu\|x\|_2\f$.
     */
    Norm2();

    /**
     * Constructs an instance of the function \f$f(x) = \mu\|x\|_2\f$
     * @param mu
     */
    explicit Norm2(double mu);

    /**
     * Default destructor.
     */
    virtual ~Norm2();

    virtual int dualNorm(Matrix& x, double& norm);

    virtual int call(Matrix& x, double& f);

    virtual int callProx(Matrix& x, double gamma, Matrix& prox);

    virtual int callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    virtual FunctionOntologicalClass category();

private:

    double m_mu;

};

#endif	/* NORM2_H */

