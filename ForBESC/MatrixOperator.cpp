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

