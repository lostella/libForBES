/* 
 * File:   Norm1.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on October 30, 2015, 4:05 PM
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

#include "Norm1.h"
#include <cmath>

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

int Norm1::callProx(Matrix& x, double gamma, Matrix& prox) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    double gm = gamma * m_mu;
    for (size_t i = 0; i < x.getNrows(); i++) {
        double xi;
        xi = x.get(i, 0);
        if (xi >= gm) {
            prox.set(i, 0, xi - gm);
        } else if (xi <= -gm) {
            prox.set(i, 0, xi + gm);
        }
    }
    return ForBESUtils::STATUS_OK;
}

int Norm1::callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    double gm = gamma * m_mu;
    f_at_prox = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        double pi;
        double xi = x.get(i, 0);
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

int Norm1::dualNorm(Matrix& x, double& norm) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    norm = std::abs(x.get(0, 0));
    for (size_t i = 1; i < x.getNrows(); i++) {
        double absi = std::abs(x.get(i, 0));
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





