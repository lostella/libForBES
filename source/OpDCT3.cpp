/* 
 * File:   OpDCT3.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 15, 2015, 6:40 PM
 */

#include "OpDCT3.h"

OpDCT3::OpDCT3(size_t m_dimension) :
LinearOperator(), m_dimension(m_dimension) {
}

OpDCT3::OpDCT3() {
    m_dimension = 0;
}

OpDCT3::~OpDCT3() {
}

Matrix OpDCT3::call(Matrix& x) {
    size_t n = x.length();
    Matrix Tx(n, 1);
    double yk;
    double x0_2 = x.get(0, 0) / 2.0;
    for (size_t k = 0; k < n; k++) {
        yk = x0_2;
        for (size_t i = 1; i < n; i++) {
            yk += (x.get(i, 0) * std::cos(i * M_PI * (k + 0.5) / n));
        }
        Tx.set(k, 0, yk);
    }
    return Tx;
}

Matrix OpDCT3::callAdjoint(Matrix& x) {
    Matrix s;
    return s;
}

size_t OpDCT3::dimensionIn() {
    return 0;
}

size_t OpDCT3::dimensionOut() {
    return 0;
}

bool OpDCT3::isSelfAdjoint() {
    return false;
}






