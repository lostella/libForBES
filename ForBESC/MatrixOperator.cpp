/* 
 * File:   MatrixOperator.cpp
 * Author: chung
 * 
 * Created on July 24, 2015, 7:31 PM
 */

#include "MatrixOperator.h"

Matrix& MatrixOperator::GetMatrix() const {
    return A;
}

bool MatrixOperator::isSelfAdjoint() {
    return m_isSelfAdjoint;
}

void MatrixOperator::SetMatrix(Matrix& A) {
    this->A = A;
    m_isSelfAdjoint = A.isSymmetric();
}

MatrixOperator::MatrixOperator(Matrix& A) :
A(A) {
    if (A.isSymmetric()) {
        this->m_isSelfAdjoint = true;
    }
}

MatrixOperator::~MatrixOperator() {
}

Matrix MatrixOperator::call(Matrix& x) {
    Matrix y;
    y = A*x;
    return y;
}

Matrix MatrixOperator::callAdjoint(Matrix& x) {
    if (isSelfAdjoint()) {
        return call(x);
    }
    Matrix y;
    A.transpose();
    y = A*x;
    A.transpose();
    return y;
}

size_t MatrixOperator::dimensionIn() {
    return A.getNcols();
}

size_t MatrixOperator::dimensionOut() {
    return A.getNrows();
}


