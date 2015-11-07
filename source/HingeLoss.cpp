/* 
 * File:   HingeLoss.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on October 29, 2015, 10:49 PM
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

#include "HingeLoss.h"
#include <iostream>

HingeLoss::HingeLoss(Matrix& b, double mu) :
Function(), m_mu(mu) {
    m_b = &b;
}

HingeLoss::HingeLoss(Matrix& b) :
Function() {
    m_b = &b;
    m_mu = 1.0;
}

HingeLoss::~HingeLoss() {
}

int HingeLoss::call(Matrix& x, double& f) {
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    f = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        double si = 1 - m_b->get(i, 0) * x.get(i, 0);
        if (si > 0) {
            f += si;
        }
    }
    f *= m_mu;
    return ForBESUtils::STATUS_OK;
}

int HingeLoss::callProx(Matrix& x, double gamma, Matrix& prox) {
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    double gm = gamma*m_mu;
    for (size_t i = 0; i < x.getNrows(); i++) {
        double bi;
        double bxi;
        bi = m_b->get(i, 0);
        bxi = bi * x.get(i, 0);
        if (bxi < 1) {
            prox.set(i, 0, bi * std::min(1.0, bxi + gm));
        } else {
            prox.set(i, 0, x.get(i, 0));
        }
    }
    return ForBESUtils::STATUS_OK;
}

int HingeLoss::callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    double gm = gamma*m_mu;
    f_at_prox = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        double si;
        double pi;
        double bi;
        double bxi;
        bi = m_b->get(i, 0);
        bxi = bi * x.get(i, 0);
        if (bxi < 1) {
            pi = bi * std::min(1.0, bxi + gm);
        } else {
            pi = x.get(i, 0);
        }
        si = 1 - m_b->get(i, 0) * pi;
        if (si > 0) {
            f_at_prox += si;
        }
        prox.set(i, 0, pi);
    }
    f_at_prox *= m_mu;
    return ForBESUtils::STATUS_OK;
}

FunctionOntologicalClass HingeLoss::category() {
    FunctionOntologicalClass hingeLoss("HingeLoss");
    hingeLoss.set_defines_f(true);
    hingeLoss.set_defines_grad(true);
    hingeLoss.set_defines_prox(true);
    hingeLoss.getSuperclasses().push_back(FunctionOntologyRegistry::loss());
    return hingeLoss;
}
