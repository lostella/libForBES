/* 
 * File:   DistanceToBall2.cpp
 * Author: chung
 * 
 * Created on January 14, 2016, 2:02 AM
 */

#include <complex>

#include "DistanceToBall2.h"

DistanceToBall2::DistanceToBall2() {
    m_center = NULL;
    m_rho = 1;
    m_w = 1;
}

DistanceToBall2::DistanceToBall2(double w)
: Function(), m_w(w) {
    m_rho = 1;
    m_center = NULL;
}

DistanceToBall2::DistanceToBall2(double w, double rho)
: Function(), m_rho(rho), m_w(w) {
    m_center = NULL;
}

DistanceToBall2::DistanceToBall2(double w, double rho, Matrix& center)
: Function(), m_rho(rho), m_w(w), m_center(&center) {

}

DistanceToBall2::~DistanceToBall2() {
}

int DistanceToBall2::call(Matrix& x, double& f, Matrix& grad) {
    double norm_x_minus_c;
    f = 0.0;
    grad = x;
    if (m_center != NULL) {
        grad -= *m_center;
    }
    norm_x_minus_c = std::sqrt((grad * grad)[0]);
    if (norm_x_minus_c > m_rho) {
        f = 0.5 * m_w * std::pow(norm_x_minus_c - m_rho, 2);
        grad *= (1-m_rho/norm_x_minus_c);
    } else {
        for (size_t i = 0; i < grad.getNrows(); i++){
            grad[i] = 0.0;
        }
    }    
    if (std::abs(m_w - 1.0) > 1e-14) {
        grad *= m_w;
    }
    return ForBESUtils::STATUS_OK;
}

int DistanceToBall2::call(Matrix& x, double& f) {
    double norm_x_minus_c = 0.0;
    f = 0.0;
    if (m_center != NULL) {
        for (size_t i = 0; i < x.getNrows(); i++) {
            norm_x_minus_c += std::pow(x[i] - m_center->get(i), 2);
        }
        norm_x_minus_c = std::sqrt(norm_x_minus_c);
    } else {
        norm_x_minus_c = std::sqrt((x * x)[0]);
    }
    if (norm_x_minus_c > m_rho) {
        f = 0.5 * m_w * std::pow(norm_x_minus_c - m_rho, 2);
    }
    return ForBESUtils::STATUS_OK;
}

FunctionOntologicalClass DistanceToBall2::category() {
    return FunctionOntologyRegistry::indicator();
}





