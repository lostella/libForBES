/* 
 * File:   MatrixOperator.cpp
 * Author: chung
 * 
 * Created on July 24, 2015, 7:31 PM
 */

#include "MatrixOperator.h"

MatrixOperator::~MatrixOperator() {
}

Matrix MatrixOperator::call(Matrix& x) {
    Matrix y;
    y = A*x;
    return y;
}

size_t MatrixOperator::dimensionIn() {
    return A.getNcols();
}

size_t MatrixOperator::dimensionOut() {
    return A.getNrows();
}


