/* 
 * File:   IndBall2.h
 * Author: Pantelis Sopasakis
 *
 * Created on November 3, 2015, 4:18 PM
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

#ifndef INDBALL2_H
#define	INDBALL2_H

#include "Function.h"


/**
 * \class IndBall2
 * \brief Indicator of a Euclidean ball
 * \version version 0.1
 * \ingroup Functions
 * \date Created on November 3, 2015, 4:18 PM
 * \author Pantelis Sopasakis
 * 
 * This class extends the class Function and implements the indicator function of
 * a Euclidean ball, that is a set of the form
 * 
 * \f[
 *  \mathcal{B}_{\rho}(x_c) = \{x \mid \|x-x_c\|_2 \leq \rho \},
 * \f]
 * 
 * for some vector \f$x_c\in\mathbb{R}^n\f$ and scalar \f$\rho>0\f$. The indicator
 * of the Euclidean ball is the function
 * 
 * \f[
 *  g(x) = \delta(x \mid \mathcal{B}_{\rho}(x_c)) = 
 * \begin{cases}
 * 0,&\text{if } \|x-x_c\|\leq \rho\\
 * \infty,&\text{otherwise}
 * \end{cases}
 * \f]
 * 
 * The proximal operator of this function is computed as follows
 * 
 * \f[
 * \mathrm{prox}_{\gamma g}(v) = \mathrm{proj}(v|\mathcal{B}_{\rho}(x_c)) = \begin{cases}
 * v,&\text{if } \|x-x_c\|\leq \rho\\
 * x_c+\rho\frac{x-x_c}{\|x-x_c\|},&\text{otherwise}
 * \end{cases}
 * \f]
 */
class IndBall2 : public Function {
public:
    /**
     * Construct a new instance of IndBall2 centered at the origin \f$x_c=0\f$ and
     * with radius \f$\rho=1.0\f$.
     */
    IndBall2();
    
    /**
     * Construct an instance of IndBall2 centered at the origin and with a given
     * radius
     * @param rho the radius of the ball
     */
    explicit IndBall2(double rho);
    
    /**
     * Create an instance of IndBall2 for a Euclidean ball centered at a given
     * point and with a given radius.
     * @param rho radius of the ball
     * @param c Center of the ball
     */
    IndBall2(double rho, Matrix& c);
    
    virtual ~IndBall2();    

    virtual int callProx(const Matrix& x, double gamma, Matrix& prox);

    virtual int callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    virtual FunctionOntologicalClass category();

    
private:
    Matrix * m_xc;
    double m_rho;
    bool m_is_xc_zero;

    double norm_div(const Matrix& x);
};

#endif	/* INDBALL2_H */

