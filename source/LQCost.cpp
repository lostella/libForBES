/* 
 * File:   LQCost.cpp
 * Author: chung
 * 
 * Created on March 3, 2016, 2:00 PM
 */

#include "LQCost.h"

LQCost::LQCost() {
    nullify_all();
}

int LQCost::factor_step() {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

void LQCost::nullify_all() {
    m_A = NULL;
    m_B = NULL;
    m_f = NULL;
    m_Q = NULL;
    m_QN = NULL;
    m_R = NULL;
    m_S = NULL;
    m_r = NULL;
    m_q = NULL;
    m_qN = NULL;
    m_N = 0;
    m_p = NULL;
    m_L = NULL;
    m_K = NULL;
    m_Rbar = NULL;
    m_d = NULL;
    m_s = NULL;
    m_RbarFactor = NULL;
}

LQCost::~LQCost() {
    if (m_L != NULL) {
        delete m_L;
        m_L = NULL;
    }
    if (m_K != NULL) {
        delete m_K;
        m_K = NULL;
    }
    if (m_Rbar != NULL) {
        delete m_Rbar;
        m_Rbar = NULL;
    }
    if (m_d != NULL) {
        delete m_d;
        m_d = NULL;
    }
    if (m_s != NULL) {
        delete m_s;
        m_s = NULL;
    }
    if (m_RbarFactor != NULL) {
        delete m_RbarFactor;
        m_RbarFactor = NULL;
    }
}

FunctionOntologicalClass LQCost::category() {
    FunctionOntologicalClass meta("QuadraticOverAffine");
    meta.set_defines_conjugate(true);
    meta.set_defines_conjugate_grad(true);
    meta.add_superclass(FunctionOntologyRegistry::conj_quadratic());
    return meta;
}

int LQCost::callConj(Matrix& x, double& f_star, Matrix& grad) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

int LQCost::callConj(Matrix& x, double& f_star) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}


