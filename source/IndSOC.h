/*
 * File:   IndSOC.h
 * Author: Lorenzo Stella
 *
 * Created on September 21, 2015, 10:07 AM
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

#ifndef INDSOC_H
#define	INDSOC_H

#include "Matrix.h"
#include "Function.h"
#include <algorithm>
#include <math.h>

/**
 * \class IndSOC
 * \brief %Indicator of a second-order cone
 * \version 0.0
 * \author Lorenzo Stella
 * \date Created on September 21, 2015, 10:07 AM
 * 
 * This class implements the indicator function of second order cones (SOC).
 * If \f$n\f$ is the dimension of the considered space, then this is the set
 *
 * \f[
 *    \mathrm{SOC}(n) = \left\{(x, t) \in \mathbb{R}^n : 
 *      x \in \mathbb{R}^{n-1} \text{and } \|x\| \leq  t\right\}
 * \f]
 *
 * The dimension \f$n\f$ must be provided to the constructor.
 * 
 * \ingroup Functions
 */
class IndSOC : public Function {
public:
    
    using Function::call;

    /**
     * Constructor for instances of IndSOC given the dimension \f$n\f$
     * of the cone.
     * 
     * @param n dimension of the space in which the second-order cone is
     * embedded. 
     */
    explicit IndSOC(int n);

    virtual ~IndSOC();

    virtual int call(Matrix& x, double& f);

    virtual int callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    virtual int callProx(Matrix& x, double gamma, Matrix& prox);
    
    virtual FunctionOntologicalClass category();


protected:

private:

    int m_n;

};

#endif	/* INDSOC_H */

