/* 
 * File:   QuadraticOperator.cpp
 * Author: chung
 * 
 * Created on July 24, 2015, 8:49 PM
 */

#include "QuadraticOperator.h"
#include "ForBESUtils.h"

QuadraticOperator::~QuadraticOperator() {
}

int QuadraticOperator::call(Matrix& x, double& f) {
    Matrix Tx = T.call(x);  /* Tx = T[x]; assertion: Tx is a vector of the same length as x */    
    if (Tx.isEmpty()){
        throw std::logic_error("Linear operator returned an empty vector");
    }
    if (!Tx.isColumnVector()){
        throw std::logic_error("Linear operator did not return a column vector");
    }
    if (Tx.getNrows() != x.getNrows()){
        throw std::logic_error("Incompatible dimensions between T(x) and x");
    }    
    f = (x * Tx).get(0, 0);
    return ForBESUtils::STATUS_OK;
}

int QuadraticOperator::category() {
    return Function::CAT_QUADRATIC;
}
