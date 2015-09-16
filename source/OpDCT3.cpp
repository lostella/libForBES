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
    if (m_dimension != 0 && n != m_dimension) {
        throw std::invalid_argument("x-dimension is invalid");
    }
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

Matrix OpDCT3::callAdjoint(Matrix& y) {
    size_t n = y.length();
    if (m_dimension != 0 && n != m_dimension) {
        throw std::invalid_argument("x-dimension is invalid");
    }
    Matrix Tstar_y(n, 1);
    double tk = 0.0;
    for (size_t i = 0; i < n; i++) {
        tk += y.get(i, 0) / 2.0;
    }
    Tstar_y.set(0, 0, tk);
    for (size_t k = 1; k < n; k++) {
        tk = 0.0;
        for (size_t i = 0; i < n; i++) {
            tk += (y.get(i, 0) * std::cos(k * M_PI * (i + 0.5) / n));
        }
        Tstar_y.set(k, 0, tk);
    }

    return Tstar_y;
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






