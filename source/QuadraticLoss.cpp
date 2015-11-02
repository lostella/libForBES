/* 
 * File:   QuadraticLoss.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on October 29, 2015, 5:47 PM
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

#include "QuadraticLoss.h"

QuadraticLoss::QuadraticLoss() {
    m_is_uniform_weights = true;
    m_is_zero_p = true;
    m_uniform_w = 1.0;
    m_w = NULL;
    m_p = NULL;
}

QuadraticLoss::QuadraticLoss(double w) :
Function(), m_uniform_w(w) {
    m_is_uniform_weights = true;
    m_is_zero_p = true;
    m_w = NULL;
    m_p = NULL;
}

QuadraticLoss::QuadraticLoss(Matrix* w, Matrix* p) :
Function(), m_p(p), m_w(w) {
    if (w == NULL || p == NULL) {
        throw std::invalid_argument("Arguments cannot be NULL in this constructor");
    }
    if (!w->isColumnVector() || !p->isColumnVector()) {
        throw std::invalid_argument("Arguments w and p must be column-vectors");
    }
    if (w->getNrows() != p->getNrows()) {
        throw std::invalid_argument("w and p must be of equal size");
    }
    m_is_uniform_weights = false;
    m_is_zero_p = false;
}

QuadraticLoss::~QuadraticLoss() {
}

int QuadraticLoss::call(Matrix& x, double& f) {
    f = 0.0;
    double fi = 0.0;
    for (size_t j = 0; j < x.getNrows(); j++) {
        fi = x.get(j, 0);
        if (!m_is_zero_p) {
            fi -= m_p->get(j, 0);
        }
        fi *= fi;
        if (!m_is_uniform_weights) {
            fi *= m_w->get(j, 0);
        }
        f += fi;
    }
    if (m_is_uniform_weights) {
        f *= m_uniform_w;
    }
    f /= 2.0;
    return ForBESUtils::STATUS_OK;
}

int QuadraticLoss::callConj(const Matrix& x, double& f_star) {
    f_star = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        f_star += x.get(i, 0)*(2.0 * (m_is_zero_p ? 0.0 : m_p->get(i, 0))
                + (x.get(i, 0) / (m_is_uniform_weights ? m_uniform_w : m_w->get(i, 0))));
    }
    f_star /= 2.0;
    return ForBESUtils::STATUS_OK;
}

int QuadraticLoss::callConj(const Matrix& x, double& f_star, Matrix& grad) {
    double gradi = 0.0;
    double pi = 0.0;
    f_star = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        pi = m_is_zero_p ? 0.0 : m_p->get(i, 0);
        gradi = pi + x.get(i, 0) / (m_is_uniform_weights ? m_uniform_w : m_w->get(i, 0));
        grad.set(i, 0, gradi);
        f_star += x.get(i, 0)*(gradi + pi);
    }
    f_star /= 2.0;
    return ForBESUtils::STATUS_OK;
}

FunctionOntologicalClass QuadraticLoss::category() {
    FunctionOntologicalClass quadLoss("QuadraticLoss");
    quadLoss.set_defines_f(true);
    quadLoss.set_defines_conjugate(true);
    quadLoss.set_defines_conjugate_grad(true);
    quadLoss.set_defines_grad(false);
    quadLoss.getSuperclasses().push_back(FunctionOntologyRegistry::loss());
    return quadLoss;
}




