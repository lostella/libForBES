/* 
 * File:   IndPos.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on July 26, 2015, 4:36 PM
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

#include <math.h>
#include <cmath>

#include "IndPos.h"

IndPos::~IndPos() {
}

IndPos::IndPos(double& uniform_lb) : m_uniform_lb(&uniform_lb) {
    m_lb = NULL;
}

IndPos::IndPos(Matrix& lb) : m_lb(&lb) {
    m_uniform_lb = NULL;
}

IndPos::IndPos() {
    m_uniform_lb = NULL;
    m_lb = NULL;
}

int IndPos::call(Matrix& x, double& f) {
    f = 0.0;
    if (m_uniform_lb != NULL) {
        for (size_t i = 0; i < x.getNrows(); i++) {
            if (x[i] < *m_uniform_lb) {
                f = INFINITY;
                break;
            }
        }
    } else if (m_lb != NULL) {
        for (size_t i = 0; i < x.getNrows(); i++) {
            if (x[i] < m_lb->get(i)) {
                f = INFINITY;
                break;
            }
        }
    } else {
        for (size_t i = 0; i < x.getNrows(); i++) {
            if (x[i] < 0.0) {
                f = INFINITY;
                break;
            }
        }
    }
    return ForBESUtils::STATUS_OK;
}

int IndPos::callConj(Matrix& y, double& f_star) {
    f_star = 0.0;
    if (m_uniform_lb != NULL) {
        for (size_t i = 0; i < y.getNrows(); i++) {
            if (y[i] > 0) {
                f_star = INFINITY;
                break;
            }
            f_star += y[i];
        }
        if (!std::isinf(f_star)) {
            f_star *= *m_uniform_lb;
        }
    } else if (m_lb != NULL) {
        for (size_t i = 0; i < y.getNrows(); i++) {
            if (y[i] > 0) {
                f_star = INFINITY;
                break;
            }
            f_star += y[i] * m_lb->get(i);
        }
    } else { /* this is like m_uniform_lb = 0.0 */
        for (size_t i = 0; i < y.getNrows(); i++) {
            if (y[i] > 0) {
                f_star = INFINITY;
                break;
            }
        }
    }


    return ForBESUtils::STATUS_OK;
}

int IndPos::callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    f_at_prox = 0.0;
    if (m_uniform_lb != NULL) {
        for (size_t i = 0; i < x.getNrows(); i++) {
            prox[i] = std::max(x[i], *m_uniform_lb);
        }
    } else if (m_lb != NULL) {
        for (size_t i = 0; i < x.getNrows(); i++) {
            prox[i] = std::max(x[i], m_lb->get(i));
        }
    } else {
        x.plusop(&prox);
    }
    return ForBESUtils::STATUS_OK;
}

int IndPos::callProx(Matrix& x, double gamma, Matrix& prox) {
    double val;
    return callProx(x, gamma, prox, val);
}

FunctionOntologicalClass IndPos::category() {
    FunctionOntologicalClass cat = FunctionOntologyRegistry::indicator();
    cat.set_defines_conjugate(true);
    return cat;
}




