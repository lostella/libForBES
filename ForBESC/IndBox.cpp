/* 
 * File:   IndBox.cpp
 * Author: chung
 * 
 * Created on July 26, 2015, 5:22 PM
 */

#include "IndBox.h"

IndBox::IndBox(double& uniform_lb, double& uniform_ub) : Function() {
    this->m_uniform_lb = &uniform_lb;
    this->m_uniform_ub = &uniform_ub;
}

IndBox::IndBox(Matrix& lb, Matrix& ub) : Function() {
    /* empty */
    if (!lb.isColumnVector()) {
        throw std::invalid_argument("LB must be a vector");
    }
    if (!ub.isColumnVector()) {
        throw std::invalid_argument("UB must be a vector");
    }
    m_lb = &lb;
    m_ub = &ub;
}

IndBox::~IndBox() {
}

int IndBox::category() {
    return Function::CAT_INDICATOR;
}

int IndBox::call(Matrix& x, double& f) {
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a vector");
    }
    bool isInside = true;
    size_t i = 0;

    if (m_uniform_lb != NULL && !isinf(-(*m_uniform_lb))) { /* there's a uniform LB and this is not -inf*/
        while (i < x.getNrows() && isInside) {
            if (m_uniform_lb != NULL) {
                isInside = isInside && (x.get(i, 0) >= *m_uniform_lb);
            }
            i++;
        }
        i = 0;
    }

    if (m_uniform_ub != NULL && !isinf(*m_uniform_ub)) {/* there's a uniform UB and this is not +inf*/
        while (i < x.getNrows() && isInside) {
            if (m_uniform_lb != NULL) {
                isInside = isInside && x.get(i, 0) <= *m_uniform_ub;
            }
            i++;
        }
        i = 0;
    }

    if (m_lb != NULL) {
        while (i < x.getNrows() && isInside) {
            isInside = isInside && x.get(i, 0) >= m_lb->get(i, 0);
            i++;
        }
        i = 0;
    }

    if (m_ub != NULL) {
        while (i < x.getNrows() && isInside) {
            isInside = isInside && x.get(i, 0) <= m_ub->get(i, 0);
            i++;
        }
    }
    f = isInside ? 1.0 : INFINITY;
    return ForBESUtils::STATUS_OK;
}

/* PROTECTED METHODS */

void IndBox::SetLb(Matrix* lb) {
    m_lb = lb;
}

void IndBox::SetUb(Matrix* ub) {
    m_ub = ub;
}

void IndBox::SetUniform_lb(double* uniform_lb) {
    m_uniform_lb = uniform_lb;
}

void IndBox::SetUniform_ub(double* uniform_ub) {
    m_uniform_ub = uniform_ub;
}

Matrix* IndBox::GetLb() const {
    return m_lb;
}

Matrix* IndBox::GetUb() const {
    return m_ub;
}

double* IndBox::GetUniform_lb() const {
    return m_uniform_lb;
}

double* IndBox::GetUniform_ub() const {
    return m_uniform_ub;
}
