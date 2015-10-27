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
 * \brief A linear operator T(x) = M*x, where M is a \c Matrix
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on July 24, 2015, 7:31 PM
 * \ingroup LinOp
 */
class MatrixOperator : public LinearOperator {
public:

    MatrixOperator(Matrix& A);
    
    Matrix& GetMatrix() const;

    void SetMatrix(Matrix& A);

    virtual bool isSelfAdjoint();

    virtual Matrix call(Matrix& x);

    virtual Matrix callAdjoint(Matrix& x);

    virtual size_t dimensionIn();

    virtual size_t dimensionOut();

    virtual ~MatrixOperator();

private:
    Matrix &A;
    bool m_isSelfAdjoint;
};

#endif	/* MATRIXOPERATOR_H */

