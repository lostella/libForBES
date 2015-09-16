/* 
 * File:   OpGradient.cpp
 * Author: chung
 * 
 * Created on September 16, 2015, 1:35 AM
 */

#include "OpGradient.h"

OpGradient::OpGradient() {
    m_dimension = 0;
}

OpGradient::OpGradient(size_t n) : LinearOperator(), m_dimension(n) {
}

OpGradient::~OpGradient() {
}

Matrix OpGradient::call(Matrix& x) {
    const size_t n = x.length();
    if (m_dimension != 0 && n != m_dimension) {
        throw std::invalid_argument("x-dimension is invalid");
    }
    Matrix Tx(n - 1, 1);
    if (n <= 1) {
        return Tx;
    }
    for (size_t i = 0; i < n - 1; i++) {
        Tx.set(i, 0, x.get(i + 1, 0) - x.get(i, 0));
    }
    return Tx;
}

Matrix OpGradient::callAdjoint(Matrix& y) {
    const size_t n = y.length() + 1;
    if (m_dimension != 0 && n != m_dimension) {
        throw std::invalid_argument("x-dimension is invalid");
    }
    Matrix Tstar_x(n, 1);
    Tstar_x.set(0, 0, -y.get(0, 0));
    for (size_t i = 1; i < n - 1; i++) {
        Tstar_x.set(i, 0, y.get(i - 1, 0) - y.get(i, 0));
    }
    Tstar_x.set(n - 1, 0, y.get(n - 2, 0));
    return Tstar_x;
}

size_t OpGradient::dimensionIn() {
    return m_dimension;
}

size_t OpGradient::dimensionOut() {
    return m_dimension == 0 ? 0 : (m_dimension - 1);
}

bool OpGradient::isSelfAdjoint() {
    return false;
}






