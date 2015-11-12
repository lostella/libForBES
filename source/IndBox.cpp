/* 
 * File:   IndBox.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on July 26, 2015, 5:22 PM
 * 
 * ForBES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *  
 * ForBES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with ForBES. If not, see <http://www.gnu.org/licenses/>.
 */


#include "IndBox.h"
#include <cmath>
#include <assert.h>

IndBox::IndBox(double& uniform_lb, double& uniform_ub) : Function() {
    m_uniform_lb = &uniform_lb;
    m_uniform_ub = &uniform_ub;
    m_lb = NULL;
    m_ub = NULL;
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
    m_uniform_lb = NULL;
    m_uniform_ub = NULL;
}

IndBox::~IndBox() {
}


int IndBox::call(Matrix& x, double& f) {
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a vector");
    }
    bool isInside = true;
    volatile size_t i = 0;

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
    f = isInside ? 0.0 : INFINITY;
    return ForBESUtils::STATUS_OK;
}

int IndBox::callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    f_at_prox = 0.0;
    assert(x.isColumnVector());

    // Step 1
    // prox = max(LB, x)
    for (size_t i = 0; i < x.getNrows(); i++) {
        if (m_lb != NULL) {
            prox.set(i, 0, std::max(m_lb->get(i, 0), x.get(i, 0)));
        } else if (m_uniform_lb != NULL) {
            prox.set(i, 0, std::max(*m_uniform_lb, x.get(i, 0)));
        }
    }

    // Step 2
    // prox = min(UB, prox)
    for (size_t i = 0; i < x.getNrows(); i++) {
        if (m_ub != NULL) {
            prox.set(i, 0, std::min(m_ub->get(i, 0), prox.get(i, 0)));
        } else if (m_uniform_ub != NULL) {
            prox.set(i, 0, std::min(*m_uniform_ub, prox.get(i, 0)));
        }
    }
    return ForBESUtils::STATUS_OK;
}

int IndBox::callProx(Matrix& x, double gamma, Matrix& prox) {
    double val;
    return callProx(x, gamma, prox, val);
}

FunctionOntologicalClass IndBox::category() {
    return FunctionOntologyRegistry::indicator();
}


/* PROTECTED METHODS */
//LCOV_EXCL_START
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
//LCOV_EXCL_STOP