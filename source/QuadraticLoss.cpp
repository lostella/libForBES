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

QuadraticLoss::QuadraticLoss(Matrix& w, Matrix& p) :
Function() {
    if (!w.isColumnVector() || !p.isColumnVector()) {
        throw std::invalid_argument("Arguments w and p must be column-vectors");
    }
    if (w.getNrows() != p.getNrows()) {
        throw std::invalid_argument("w and p must be of equal size");
    }
    m_is_uniform_weights = false;
    m_is_zero_p = false;
    m_w = &w;
    m_p = &p;
    m_uniform_w = 0.0;
}

QuadraticLoss::~QuadraticLoss() {
}

int QuadraticLoss::call(Matrix& x, double& f) {
    f = 0.0;
    for (size_t j = 0; j < x.getNrows(); j++) {
        double fi;
        fi = x[j];
        if (!m_is_zero_p) {
            fi -= m_p->getData()[j];
        }
        fi *= fi;
        if (!m_is_uniform_weights) {
            fi *= m_w->getData()[j];
        }
        f += fi;
    }
    if (m_is_uniform_weights) {
        f *= m_uniform_w;
    }
    f /= 2.0;
    return ForBESUtils::STATUS_OK;
}

int QuadraticLoss::call(Matrix& x, double& f, Matrix& grad) {
    f = 0.0;
    for (size_t j = 0; j < x.getNrows(); j++) {
        double fi;
        double gi;
        fi = x[j];
        if (!m_is_zero_p) {
            fi -= m_p->getData()[j];
        }
        gi = fi;
        fi *= fi;
        if (!m_is_uniform_weights) {
            double w = m_w->getData()[j];
            fi *= w;
            gi *= w;
        }
        f += fi;
        grad[j] = gi;
    }
    if (m_is_uniform_weights) {
        f *= m_uniform_w;
    }
    f /= 2.0;
    return ForBESUtils::STATUS_OK;
}

int QuadraticLoss::callConj(Matrix& x, double& f_star) {
    f_star = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        f_star += x[i]*(2.0 * (m_is_zero_p ? 0.0 : m_p->get(i))
                + (x[i] / (m_is_uniform_weights ? m_uniform_w : m_w->get(i))));
    }
    f_star /= 2.0;
    return ForBESUtils::STATUS_OK;
}

int QuadraticLoss::callConj(Matrix& x, double& f_star, Matrix& grad) {
    f_star = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        double gradi;
        double pi;
        pi = m_is_zero_p ? 0.0 : m_p->get(i);
        gradi = pi + x[i] / (m_is_uniform_weights ? m_uniform_w : m_w->get(i));
        grad.set(i, 0, gradi);
        f_star += x[i]*(gradi + pi);
    }
    f_star /= 2.0;
    return ForBESUtils::STATUS_OK;
}

int QuadraticLoss::hessianProduct(Matrix& x, Matrix& z, Matrix& Hz) {
    if (m_is_uniform_weights) {
        Hz = m_uniform_w * z;
    } else {
        m_w->toggle_diagonal();
        Hz = (*m_w)*z;
        m_w->toggle_diagonal();
    }
    return ForBESUtils::STATUS_OK;
}

FunctionOntologicalClass QuadraticLoss::category() {
    FunctionOntologicalClass quadLoss("QuadraticLoss");
    quadLoss.set_defines_f(true);
    quadLoss.set_defines_conjugate(true);
    quadLoss.set_defines_conjugate_grad(true);
    quadLoss.set_defines_grad(true);
    quadLoss.add_superclass(FunctionOntologyRegistry::loss());
    quadLoss.add_superclass(FunctionOntologyRegistry::quadratic());
    return quadLoss;
}
