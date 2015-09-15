/* 
 * File:   LinearOperator.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 24, 2015, 5:05 PM
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

#ifndef LINEAROPERATOR_H
#define	LINEAROPERATOR_H

#include "Matrix.h"

/**
 * \class LinearOperator
 * \brief A linear operator T(x)
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on July 24, 2015, 5:05 PM
 * 
 * An interface for an arbitrary linear operator.
 * This is an abstract class which cannot be instantiated (it involves
 * pure virtual members).
 */
class LinearOperator {
public:

    /**
     * Computes the value of the operator at a vector <code>x</code>.
     * @param x vector where the operator should be calculated
     * @return Value <code>T(x)</code>
     */
    virtual Matrix call(Matrix& x) = 0;
    
    /**
     * Calls the adjoint of this operator. 
     * 
     * For an operator \f$T:\mathbb{R}^n \to \mathbb{R}^m\f$, its adjoint is defined
     * as an operator \f$T^*:\mathbb{R}^m \to \mathbb{R}^n\f$ so that for all \f$x\in\mathbb{R}^n\f$
     * and every \f$y\in\mathbb{R}^m\f$ it is
     * \f[
     *  \langle y, T(x) \rangle = \langle T^*(y), x \rangle
     * \f]
     * For a matrix \f$M\in\mathbb{R}^{m\times n}\f$ its adjoint is its transpose,
     * that is \f$M^* = M'\f$ and, indeed, for all \f$x\in\mathbb{R}^n\f$ 
     * and \f$y\in\mathbb{R}^m\f$ it is \f$\langle y, Mx \rangle = y'Mx = (M'y)'x = \langle M'y, x \rangle = 
     * \langle M^*(y), x \rangle\f$.    
     * 
     * @param x Vector <code>x</code>
     * @return The result \f$T^*(x)\f$
     */
    virtual Matrix callAdjoint(Matrix& x) = 0;
    
    /**
     * Whether this operator is self-adjoint.
     * 
     * An operator is called <em>self-adjoint</em> if it coincides with its
     * adjoint.
     * 
     * @return <code>true</code> if and only if the operator is self-adjoint.
     */
    virtual bool isSelfAdjoint() = 0;
    
    /**
     * For a linear operator \f$T:\mathbb{R}^n \to \mathbb{R}^p\f$, it returns
     * the dimension of its domain, i.e., <code>n</code>. 
     * 
     * In some cases, this
     * dimension may be undefined, so then this method returns <code>0</code>; such a 
     * case may be a linear operator which reverses the order of the entries
     * of <code>x</code> which may be applied to vectors of any size. It is, however, 
     * strongly recommended to define the operator's dimensions in all cases.
     * 
     * @return Dimension of <code>x</code>
     */
    virtual size_t dimensionIn() = 0;
    
    /**
     * For a linear operator \f$T:\mathbb{R}^n \to \mathbb{R}^p\f$, it returns
     * the dimension of its range, i.e., <code>p</code>.
     * 
     * In some cases, this
     * dimension may be undefined, so then this method returns <code>0</code>; such a 
     * case may be a linear operator which reverses the order of the entries
     * of <code>x</code> which may be applied to vectors of any size. It is, however, 
     * strongly recommended to define the operator's dimensions in all cases.
     * 
     * @return Dimension of <code>T(x)</code>
     */
    virtual size_t dimensionOut() = 0;

    /**
     * Destructor for LinearOperator objects.
     */
    virtual ~LinearOperator();

protected:
    LinearOperator();
    LinearOperator(const LinearOperator& orig);

};

#endif	/* LINEAROPERATOR_H */

