/* 
 * File:   OpReverseVector.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 15, 2015, 12:57 PM
 */

#include "OpReverseVector.h"

OpReverseVector::OpReverseVector() {
    m_vectorDim = 0;
}

OpReverseVector::~OpReverseVector() {
}

OpReverseVector::OpReverseVector(size_t n) {
    this->m_vectorDim = n;
}

Matrix OpReverseVector::call(Matrix& x) {
    Matrix y(x);
    double temp;
    for (size_t i = 0; i < y.length() / 2; ++i) {
        temp = y.get(y.length() - i - 1, 0);
        y.set(y.length() - i - 1, 0, y.get(i, 0));
        y.set(i, 0, temp);
    }
    return y;
}

Matrix OpReverseVector::callAdjoint(Matrix& x) {
    return call(x);
}

size_t OpReverseVector::dimensionIn() {
    return m_vectorDim;
}

size_t OpReverseVector::dimensionOut() {
    return m_vectorDim;
}

bool OpReverseVector::isSelfAdjoint() {
    return true;
}







