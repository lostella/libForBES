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

/**
 * Dimension of a vector-valued operator: n x 1.
 */
#define _VECTOR_OP_DIM(N) std::pair<size_t, size_t>(static_cast<size_t> (N), static_cast<size_t> (1))

/**
 * Dimension of an empty operator: 0 x 1.
 */
#define _EMPTY_OP_DIM _VECTOR_OP_DIM(0)

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
 * 
 * Linear operators are assumed to be of the generic form \f$T:X \to Y\f$,
 * where \f$X\f$ and \f$Y\f$ are vector spaces, either \f$\mathbb{R}^n\f$ or
 * \f$\mathbb{R}^{n\times m}\f$.
 */
class LinearOperator {
public:

    /**
     * Computes the value of the operator at a vector or matrix <code>x</code>.
     * @param x vector or matrix where the operator should be calculated
     * @return Value <code>T(x)</code>
     * 
     */
    virtual Matrix call(Matrix& x) = 0;
    
    /**
     * Calls the adjoint of this operator. 
     * 
     * For an operator \f$T:X \to Y\f$, its adjoint is defined
     * as an operator \f$T^*:Y \to X\f$ so that for all \f$x\in X\f$
     * and every \f$y Y\f$ it is
     * \f[
     *  \langle y, T(x) \rangle_{Y} = \langle T^*(y), x \rangle_{X}
     * \f]
     * For a matrix \f$M\in\mathbb{R}^{m\times n}\f$ its adjoint is its transpose,
     * that is \f$M^* = M'\f$ and, indeed, for all \f$x\in\mathbb{R}^n\f$ 
     * and \f$y\in\mathbb{R}^m\f$ it is \f$\langle y, Mx \rangle = y'Mx = (M'y)'x = \langle M'y, x \rangle = 
     * \langle M^*(y), x \rangle\f$.    
     * 
     * @param x Vector or matrix <code>x</code>
     * @return The result \f$T^*(x)\f$
     */
    virtual Matrix callAdjoint(Matrix& x) = 0;
    
    /**
     * Whether this operator is self-adjoint.
     * 
     * An operator \f$T:X\to X\f$ is called <em>self-adjoint</em> if it coincides with its
     * adjoint.
     * 
     * @return <code>true</code> if and only if the operator is self-adjoint.
     */
    virtual bool isSelfAdjoint() = 0;
    
    /**
     * For a linear operator \f$T:X \to Y\f$, it returns
     * the dimension of its domain.
     * 
     * The in-dimension of an operator \f$f:\mathbb{R}^{n}\to Y\f$ is the pair <code>(n,1)</code>.
     * 
     * The in-dimension of \f$f:\mathbb{R}^{n\times m}\to Y\f$ is the pair <code>(n,m)</code>.
     * 
     * In some cases, this dimension may be undefined, so then this method returns 
     * <code>(0, ?)</code>; such a case may be a linear operator which reverses the 
     * order of the entries of <code>x</code> which may be applied to vectors of any 
     * size. In such a case, the in-dimension will be <code>(0,1)</code>.
     * It is, however, strongly recommended to define the nonzero operator's dimensions 
     * in all cases.
     * 
     * @return Dimensions of <code>x</code>
     */
    virtual std::pair<size_t, size_t> dimensionIn() = 0;
    
    /**
     * For a linear operator \f$T:X \to Y\f$, it returns
     * the dimension of its range space \f$Y\f$.
     * 
     * In some cases, this dimension may be undefined, so then this method returns 
     * <code>(0, ?)</code>; such a case may be a linear operator which reverses the 
     * order of the entries of <code>x</code> which may be applied to vectors of any 
     * size. In such a case, the in-dimension will be <code>(0,1)</code>.
     * It is, however, strongly recommended to define the nonzero operator's dimensions 
     * in all cases.
     * 
     * @return Dimensions of <code>T(x)</code>
     */
    virtual std::pair<size_t, size_t> dimensionOut() = 0;

    /**
     * Destructor for LinearOperator objects.
     */
    virtual ~LinearOperator();

protected:
    LinearOperator();

};

#endif	/* LINEAROPERATOR_H */

