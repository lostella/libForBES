/* 
 * File:   OpLinearCombination.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 14, 2015, 9:25 PM
 */

#include "OpLinearCombination.h"

OpLinearCombination::OpLinearCombination(LinearOperator& A, LinearOperator& B, double a, double b) :
LinearOperator(), A(A), B(B), a(a), b(b) {
    if (A.dimensionIn() != B.dimensionIn()) {
        throw std::invalid_argument("A and B have incompatible input dimensions");
    }
    if (A.dimensionOut() != B.dimensionOut()) {
        throw std::invalid_argument("A and B have incompatible input dimensions");
    }
}

OpLinearCombination::~OpLinearCombination() {
}

