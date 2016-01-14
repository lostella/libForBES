/* 
 * File:   MatrixOperator.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 24, 2015, 7:31 PM
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

#ifndef MATRIXOPERATOR_H
#define	MATRIXOPERATOR_H

#include "Matrix.h"
#include "LinearOperator.h"

/**
 * \class MatrixOperator
 * \brief A linear operator T(x) = M*x, where M is a %Matrix
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on July 24, 2015, 7:31 PM
 * \ingroup LinOp
 * 
 * \example matop_example.cpp
 */
class MatrixOperator : public LinearOperator {
public:

    using LinearOperator::call;
    using LinearOperator::callAdjoint;

    /**
     * Defines a constructs a new instance of MatrixOperator providing a reference
     * to an instance of Matrix.
     * @param A Matrix
     */
    explicit MatrixOperator(Matrix& A);

    /**
     * Provides access to the underlying matrix.
     * @return this operator as a matrix.
     */
    Matrix& getMatrix() const;

    /**
     * Allows the update of the underlying matrix.
     * @param A a new instance of Matrix
     */
    void setMatrix(Matrix& A);

    /**
     * Whether this operator is self-adjoint, i.e., whether the underlying matrix
     * is symmetric.
     * @return whether the operator is self-adjoint
     */
    virtual bool isSelfAdjoint();

    virtual int call(Matrix& y, double alpha, Matrix& x, double gamma);

    virtual int callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma);

    virtual std::pair<size_t, size_t> dimensionIn();

    virtual std::pair<size_t, size_t> dimensionOut();

    /**
     * Default destructor
     */
    virtual ~MatrixOperator();

private:
    Matrix & m_A; /**< matrix which defines the operator */
    bool m_isSelfAdjoint;/**< whether this is self-adjoint */
};

#endif	/* MATRIXOPERATOR_H */

