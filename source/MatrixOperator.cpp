/* 
 * File:   MatrixOperator.cpp
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

#include "MatrixOperator.h"

Matrix& MatrixOperator::GetMatrix() const {
    return m_A;
}

bool MatrixOperator::isSelfAdjoint() {
    return m_isSelfAdjoint;
}

void MatrixOperator::SetMatrix(Matrix& A) {
    this->m_A = A;
    m_isSelfAdjoint = (A.getNrows() == A.getNcols() && A.isSymmetric());
}

MatrixOperator::MatrixOperator(Matrix& A) : m_A(A) {
    if (A.isSymmetric()) {
        this->m_isSelfAdjoint = true;
    } else {
        this->m_isSelfAdjoint = false;
    }
}

MatrixOperator::~MatrixOperator() {
}

int MatrixOperator::call(Matrix& y, double alpha, Matrix& x, double gamma) {
    return Matrix::mult(y, alpha, m_A, x, gamma);
}

int MatrixOperator::callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma) {
    if (isSelfAdjoint()) {
        return call(y, alpha, x, gamma);
    }
    m_A.transpose();
    int status = Matrix::mult(y, alpha, m_A, x, gamma);    
    m_A.transpose();
    return status;
}

std::pair<size_t, size_t> MatrixOperator::dimensionIn() {
    return _VECTOR_OP_DIM(m_A.getNcols());
}

std::pair<size_t, size_t> MatrixOperator::dimensionOut() {
    return _VECTOR_OP_DIM(m_A.getNrows());
}


