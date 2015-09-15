/* 
 * File:   OpAdjoint.cpp
 * Author: chung
 * 
 * Created on September 15, 2015, 2:57 PM
 */

#include "OpAdjoint.h"

OpAdjoint::OpAdjoint(LinearOperator& op) : LinearOperator(), m_originalOperator(op) {
    
}

OpAdjoint::~OpAdjoint() {
}

Matrix OpAdjoint::call(Matrix& x) {
    return m_originalOperator.callAdjoint(x);
}

Matrix OpAdjoint::callAdjoint(Matrix& x) {
    return m_originalOperator.call(x);
}

size_t OpAdjoint::dimensionIn() {
    return m_originalOperator.dimensionOut();
}

size_t OpAdjoint::dimensionOut() {
    return m_originalOperator.dimensionIn();
}

bool OpAdjoint::isSelfAdjoint() {
    return m_originalOperator.isSelfAdjoint();
}






