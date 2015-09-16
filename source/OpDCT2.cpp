/* 
 * File:   OpDCT2.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 15, 2015, 3:36 PM
 */

#include "OpDCT2.h"

OpDCT2::OpDCT2() : LinearOperator() {
    m_dimension = 0;
}

OpDCT2::OpDCT2(size_t n) {
    m_dimension = n;
}

OpDCT2::~OpDCT2() {
}

Matrix OpDCT2::call(Matrix& x) {
    size_t n = x.length();
    Matrix Y(n, 1);
    double yk;
    for (size_t k = 0; k < n; k++) {
        yk = 0.0;
        for (size_t i = 0; i < n; i++) {
            yk += (x.get(i, 0) * std::cos(k * M_PI * (i + 0.5) / n));
        }
        Y.set(k, 0, yk);
    }
    return Y;
}

Matrix OpDCT2::callAdjoint(Matrix& x) {
    size_t n = x.length();
    Matrix Yadj(n, 1);
    double yk;
    for (size_t k = 0; k < n; k++) {
        yk = 0.0;
        for (size_t i = 0; i < n; i++) {
            yk += (x.get(i, 0) * std::cos(i * M_PI * (k + 0.5) / n));
        }
        Yadj.set(k, 0, yk);
    }
    return Yadj;
}

size_t OpDCT2::dimensionIn() {
    return m_dimension;
}

size_t OpDCT2::dimensionOut() {
    return m_dimension;
}

bool OpDCT2::isSelfAdjoint() {
    return false;
}
