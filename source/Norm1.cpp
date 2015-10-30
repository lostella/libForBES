/* 
 * File:   Norm1.cpp
 * Author: chung
 * 
 * Created on October 30, 2015, 4:05 PM
 */

#include "Norm1.h"

Norm1::Norm1() : Norm() {
    m_mu = 1.0;
}

Norm1::Norm1(double mu) : Norm(), m_mu(mu) {
    if (mu <= 0) {
        throw std::invalid_argument("Parameter mu must be positive");
    }
}

Norm1::~Norm1() {
}

int Norm1::call(Matrix& x, double& f) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    f = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        f += std::abs(x.get(i, 0));
    }
    f *= m_mu;
    return ForBESUtils::STATUS_OK;
}

int Norm1::callProx(const Matrix& x, double gamma, Matrix& prox) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    double gm = gamma * m_mu;
    double xi;
    for (size_t i = 0; i < x.getNrows(); i++) {
        xi = x.get(i, 0);
        if (xi >= gm) {
            prox.set(i, 0, xi - gm);
        } else if (xi <= -gm) {
            prox.set(i, 0, xi + gm);
        }
    }
    return ForBESUtils::STATUS_OK;
}

int Norm1::callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    double gm = gamma * m_mu;
    double xi;
    double pi;
    f_at_prox = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        xi = x.get(i, 0);
        if (xi >= gm) {
            pi = xi - gm;
        } else if (xi <= -gm) {
            pi = xi + gm;
        } else {
            pi = 0.0;
        }
        prox.set(i, 0, pi);
        f_at_prox += std::abs(pi);
    }
    f_at_prox *= m_mu;
    return ForBESUtils::STATUS_OK;
}

int Norm1::dualNorm(const Matrix& x, double& norm) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    norm = std::abs(x.get(0, 0));
    double absi;
    for (size_t i = 1; i < x.getNrows(); i++) {
        absi = std::abs(x.get(i, 0));
        if (absi > norm) {
            norm = absi;
        }
    }
    norm /= m_mu;
    return ForBESUtils::STATUS_OK;
}

FunctionOntologicalClass Norm1::category() {
    return FunctionOntologyRegistry::norm();
}





