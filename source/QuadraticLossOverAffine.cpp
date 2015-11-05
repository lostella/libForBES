/* 
 * File:   QuadraticLossOverAffine.cpp
 * Author: chung
 * 
 * Created on November 3, 2015, 3:58 PM
 */

#include "QuadraticLossOverAffine.h"
#include "LDLFactorization.h"
#include <iostream>

QuadraticLossOverAffine::QuadraticLossOverAffine(Matrix& A, Matrix& b, Matrix& w, Matrix& p) {
    m_A = &A;
    m_b = &b;
    m_w = &w;
    m_p = &p;
    if (!w.isColumnVector()) {
        throw std::invalid_argument("w is not a column vector");
    }
    /* F = A * diag(1/sqrt(w_i))_i */
    Matrix F;
    Matrix W_inv_sqrt(w);
    for (size_t i = 0; i < w.getNrows(); ++i) {
        W_inv_sqrt.set(i, 0, 1 / std::sqrt(W_inv_sqrt.get(i, 0)));
    }
    W_inv_sqrt.toggle_diagonal();
    F = A * W_inv_sqrt;
    
}

QuadraticLossOverAffine::~QuadraticLossOverAffine() {
}


