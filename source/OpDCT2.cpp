/* 
 * File:   OpDCT2.cpp
 * Author: chung
 * 
 * Created on September 15, 2015, 3:36 PM
 */

#include "OpDCT2.h"

OpDCT2::OpDCT2() : LinearOperator() {
}

OpDCT2::~OpDCT2() {
}

Matrix OpDCT2::call(Matrix& x) {
    size_t N = x.length();
    if (N == 0) {
        throw std::invalid_argument("OpDCT2::call cannot be applied to empty vectors");
    }
    Matrix Y(N, 1);
    double yk;
    for (size_t k = 0; k < N; k++) {
        yk = 0.0;
        for (size_t i = 0; i < N; i++) {
            yk += (x.get(i, 0) * std::cos(k * M_PI * (i + 0.5) / N));
        }
        Y.set(k, 0, yk);
    }
    return Y;
}

Matrix OpDCT2::callAdjoint(Matrix& x) {
    throw std::logic_error("Not implemented yet!");
}

size_t OpDCT2::dimensionIn() {
    return -1;
}

size_t OpDCT2::dimensionOut() {
    return -1;
}

bool OpDCT2::isSelfAdjoint() {
    return false;
}





