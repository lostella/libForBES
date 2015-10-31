/* 
 * File:   QuadraticLoss.h
 * Author: Pantelis Sopasakis
 *
 * Created on October 29, 2015, 5:47 PM
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

#ifndef QUADRATICLOSS_H
#define	QUADRATICLOSS_H

#include "Function.h"

/**
 * \class QuadraticLoss
 * \brief %Quadratic loss function
 * \version version 0.1
 * \ingroup Functions
 * \date Created on October 29, 2015, 5:47 PM
 * \author Pantelis Sopasakis
 * 
 * This class extends the class Function and implements the quadratic loss
 * function which is a function \f$f:\mathbb{R}^n\to\mathbb{R}\f$ defined as 
 * 
 * \f[
 *  f(x) = \frac{1}{2} \sum_{i=1}^{n} w_i(x_i - p_i)^2,
 * \f]
 *
 * where \f$w,p\in\mathbb{R}^n\f$ are given vectors.
 * 
 * Often, it is assumed that \f$w_i=w\f$ for all \f$i\f$. Then, the above is simplified as
 * 
 * \f[
 *  f(x) = \frac{w}{2} \sum_{i=1}^{n} (x_i - p_i)^2.
 * \f]
 * 
 * The conjugate of this function is differentiable and its gradient is
 * 
 * \f[
 * (\nabla f^*(y))_i = p_i + \frac{y_i}{w_i}.
 * \f]
 * 
 * while the conjugate is computed as
 * 
 * \f[
 * f^*(y) = \frac{1}{2} y' (\nabla f^*(y) + p).
 * \f]
 * 
 * 
 */
class QuadraticLoss : public Function {
public:
    /**
     * Create a new instance of QuadraticLoss assuming a uniform weight \f$w=1\f$
     * and \f$p=0\f$.
     */
    QuadraticLoss();

    /**
     * Create a new instance of QuadraticLoss using a given uniform weight \f$w\f$
     * and setting \f$p=0\f$.
     * @param w uniform weight
     */
    QuadraticLoss(double w); // p = 0

    /**
     * Create a new instance of QuadraticLoss with given parameters \f$w\f$ and \f$p\f$.
     * 
     * @param w vector of weights
     * @param p vector p
     */
    QuadraticLoss(Matrix* w, Matrix* p);

    /**
     * Destructor.
     */
    virtual ~QuadraticLoss();


    virtual int call(Matrix& x, double& f);

    virtual int callConj(const Matrix& x, double& f_star);

    virtual int callConj(const Matrix& x, double& f_star, Matrix& grad);

    virtual FunctionOntologicalClass category();




private:

    bool m_is_uniform_weights;
    bool m_is_zero_p;

    Matrix * m_w = NULL;
    Matrix * m_p = NULL;

    double m_uniform_w;

};

#endif	/* QUADRATICLOSS_H */

