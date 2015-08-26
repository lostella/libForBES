/* 
 * File:   QuadraticOperator.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on July 24, 2015, 8:49 PM
 */

#include "QuadraticOperator.h"
#include "ForBESUtils.h"

QuadraticOperator::~QuadraticOperator() {
}

int QuadraticOperator::call(Matrix& x, double& f, Matrix& grad) {
    int status = computeGradient(x, grad);
    f = (x * grad).get(0, 0)/2;
    return status;
}

int QuadraticOperator::call(Matrix& x, double& f) {
    Matrix grad;
    int status = call(x, f, grad);
    return status;
}

int QuadraticOperator::computeGradient(Matrix& x, Matrix& grad) {
    grad = T.call(x); /* Tx = T[x]; assertion: Tx is a vector of the same length as x */
    
    if (grad.isEmpty()) {
        throw std::logic_error("Linear operator returned an empty vector");
    }
    if (!grad.isColumnVector()) {
        throw std::logic_error("Linear operator did not return a column vector");
    }
    if (grad.getNrows() != x.getNrows()) {
        throw std::logic_error("Incompatible dimensions between T(x) and x");
    }
    return ForBESUtils::STATUS_OK;
}

int QuadraticOperator::category() {
    return Function::CAT_QUADRATIC;
}
