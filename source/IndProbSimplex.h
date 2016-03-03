/* 
 * File:   IndProbSimplex.h
 * Author: Pantelis Sopasakis
 *
 * Created on January 15, 2016, 2:37 AM
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

#ifndef INDPROBSIMPLEX_H
#define	INDPROBSIMPLEX_H

#include "Function.h"

/**
 * \class IndProbSimplex
 * \brief %Indicator of the probability simplex
 * \version 0.1
 * \author Pantelis Sopasakis
 * \date Created on January 15, 2016, 2:37 AM
 * 
 * This class implements the indicator of the probability simplex, that is
 * 
 * \f[
 * f(x) = \delta(x\mid P) = 
 * \begin{cases}
 *  0, &\text{if } x\in P\\
 * +\infty,&\text{otherwise}
 * \end{cases}
 * \f]
 * 
 * where \f$P\f$ is the set
 * 
 * \f[
 * P = \{x\in\mathbb{R}^n: x\geq 0, 1^{\top}x = 1\}.
 * \f]
 * 
 * The projection on this set, that is the proximal mapping of the indicator,
 * is given by
 * 
 * \f[
 * \mathrm{proj}(x\mid P) = (x-t1)_+,
 * \f]
 * 
 * where \f$(z)_+=\max \{z, 0\}\f$ (element-wise) and \f$\max_i x_i - 1 \leq t \leq \max_i x_i\f$ 
 * is a scalar such that
 * 
 * \f[
 * 1^{\top}(x-t1)_+ = 1 \Leftrightarrow \sum_{i=1}^{n}(x_i-t)_+ = 1.
 * \f]
 * 
 * This can be determined by bisection.
 */
class IndProbSimplex : public Function {
public:
    using Function::call;    
    
    /**
     * Constructs a new instance of the indicator function of the probability
     * simplex.
     * 
     */
    IndProbSimplex();

    /**
     * Default destructor.
     */
    virtual ~IndProbSimplex();
    
    virtual int call(Matrix& x, double& f);

    virtual int callProx(Matrix& x, double gamma, Matrix& prox);
    
    virtual int callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    virtual FunctionOntologicalClass category();


private:
 

};

#endif	/* INDPROBSIMPLEX_H */

