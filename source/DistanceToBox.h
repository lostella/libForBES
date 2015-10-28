/* 
 * File:   DistanceToBox.h
 * Author: Pantelis Sopasakis
 *
 * Created on October 27, 2015, 6:05 PM
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

#ifndef DISTANCETOBOX_H
#define	DISTANCETOBOX_H

#include "Function.h"

/**
 * \class DistanceToBox
 * \brief Distance from a box in R^n
 * \version version 0.1
 * \ingroup Functions
 * \date Created on October 27, 2015, 6:05 PM
 * \author Pantelis Sopasakis
 * 
 * This class extends the class Function and implements the squared Euclidean distance of
 * a vector from a box.
 * 
 * A box is a set of the form \f$B_{[l,u]}=\{x \in \mathbb{R}^n: l\leq x \leq u\}\f$
 * and the squared distance-to-box function is defined as
 * \f[
 * f(x) = \frac{1}{2}\mathrm{dist}^2 (x\mid B_{[l,u]}) = \frac{1}{2}\min_{y\in B_{[l,u]}}\|y-x\|^2_W,
 * \f]
 * where \f$W\f$ is a diagonal matrix with positive diagonal entries, and \f$\|\cdot\|_W^2\f$
 * is defined as
 * 
 * \f[
 *  \|z\|_W^2 = \sum_{i=1}^{N} W_{ii} z_i^2.
 * \f]
 * 
 * 
 * To compute this function, we use the fact that the projection of a \f$x\in\mathbb{R}^n\f$
 * on \f$B_{[l,u]}\f$ can be easily determined as
 * 
 * \f[
 *  \Pi(x\mid B_{[l,u]}) = \max (\min(x,ub), lb).
 * \f]
 * 
 * Then, we define
 * 
 * \f[
 *  d(x) = x - \Pi(x\mid B_{[l,u]}),
 * \f]
 * 
 * and we have
 * 
 * \f[
 *  \nabla f(x) = W d(x),
 * \f]
 * 
 * and
 * 
 * \f[
 *  f(x) = \frac{1}{2}\|d(x)\|_W^2.
 * \f]
 */
class DistanceToBox : public Function {
public:


    virtual ~DistanceToBox();

    /**
     * Construct a new instance of DistanceToBox with given lower and upper bounds
     * for the box and weights.
     * 
     * @param lb lower bound of the box
     * @param ub upper bound of the box
     * @param weights weights (all positive)
     */
    DistanceToBox(Matrix* lb, Matrix* ub, Matrix* weights);


    /**
     * Constructor for instances of DistanceToBox objects given the upper and lower
     * bounds for the box and a single weight (scalar, double) assuming all weights
     * are equal. In such a case, the Distance to box function can be written as
     * 
     * \f[
     *  f(x) = \frac{w}{2}\min_{y\in B_{[l,u]}}\|y-x\|^2.
     * \f]
     * 
     * @param lb lower bound of the box
     * @param ub upper bound of the box
     * @param weight scalar weight
     */
    DistanceToBox(Matrix* lb, Matrix* ub, double weight);


    virtual int call(Matrix& x, double& f, Matrix& grad);

    virtual int call(Matrix& x, double& f);

    virtual int category();




private:
    
    int compute_dx(Matrix& x, Matrix& dx) const;

    bool m_is_weights_equal; /**< Whether all weights are equal to each other. */
    bool m_is_bounds_uniform; /**< Whether box bounds are uniform. */
    double m_weight; /**< The single scalar weight if <code>m_is_weights_equal</code> is true. */
    double m_uniform_lb; /**< A uniform lower bound (with <code>m_is_bounds_uniform</code> set to true). */
    double m_uniform_ub; /**< A uniform upper bound (with <code>m_is_bounds_uniform</code> set to true). */
    Matrix* m_weights = NULL; /**< A vector of weights (if <code>m_is_weights_equal == false</code>). */
    Matrix* m_lb = NULL; /**< The box lower bound. */
    Matrix* m_ub = NULL; /**< The box upper bound. */


};

#endif	/* DISTANCETOBOX_H */

