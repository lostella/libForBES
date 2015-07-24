/* 
 * File:   MatrixOperator.h
 * Author: chung
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

class MatrixOperator : public LinearOperator {
public:

    MatrixOperator(Matrix& A) :
    A(A) {        
    }

    MatrixOperator(const MatrixOperator& other) :
    A(other.A) {
    }

    Matrix& GetMatrix() const {
        return A;
    }

    void SetMatrix(Matrix& A) {
        this->A = A;
    }


    virtual Matrix call(Matrix& x);

    virtual size_t dimensionIn();

    virtual size_t dimensionOut();

    virtual ~MatrixOperator();

private:
    Matrix &A;
};

#endif	/* MATRIXOPERATOR_H */

