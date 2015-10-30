/* 
 * File:   HuberLoss.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on October 30, 2015, 1:57 AM
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

#include "HuberLoss.h"

HuberLoss::HuberLoss(double delta) :
Function(), m_delta(delta) {
}

HuberLoss::~HuberLoss() {
}

int HuberLoss::call(Matrix& x, double& f) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    f = 0.0;
    double xi;
    for (size_t i = 0; i < x.getNrows(); i++) {
        xi = x.get(i, 0);
        if (std::abs(xi) <= m_delta) {
            f += std::pow(xi, 2) / (2.0 * m_delta);
        } else {
            f += std::abs(xi) - m_delta / 2.0;
        }
    }
    return ForBESUtils::STATUS_OK;
}

int HuberLoss::call(Matrix& x, double& f, Matrix& grad) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    f = 0.0;
    double xi;
    for (size_t i = 0; i < x.getNrows(); i++) {
        xi = x.get(i, 0);
        if (std::abs(xi) <= m_delta) {
            f += std::pow(xi, 2) / (2.0 * m_delta);
            grad.set(i, 0, xi / m_delta);
        } else {
            f += std::abs(xi) - m_delta / 2.0;
            grad.set(i, 0, std::signbit(xi) ? -1.0 : 1.0);
        }
    }
    return ForBESUtils::STATUS_OK;
}

