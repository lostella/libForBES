/* 
 * File:   OpComposition.cpp
 * Author: chung
 * 
 * Created on September 14, 2015, 9:26 PM
 */

#include "OpComposition.h"

OpComposition::~OpComposition() {
}

Matrix OpComposition::call(Matrix& x) {
    Matrix y = B.call(x);   // y = B(x)
    return A.call(y);       // z = A(y)
}

size_t OpComposition::dimensionIn() {
    return B.dimensionIn();
}

size_t OpComposition::dimensionOut() {
    return A.dimensionOut();
}

bool OpComposition::isSelfAdjoint() {
    return A.isSelfAdjoint() && B.isSelfAdjoint();
}

Matrix OpComposition::callAdjoint(Matrix& x) {
    Matrix y = A.callAdjoint(x); // y = B(x)
    return B.callAdjoint(y); // z = A(y)
}
