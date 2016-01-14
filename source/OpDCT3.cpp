/* 
 * File:   OpDCT3.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 15, 2015, 6:40 PM
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

#include "OpDCT3.h"
#include <cmath>
#include <limits>

OpDCT3::OpDCT3(size_t dimension) :
LinearOperator(),
m_dimension(_VECTOR_OP_DIM(dimension)) {
}

OpDCT3::OpDCT3() : m_dimension(_EMPTY_OP_DIM) {

}

OpDCT3::~OpDCT3() {
}

int OpDCT3::call(Matrix& y, double alpha, Matrix& x, double gamma) {
    size_t n = x.length();
    if (m_dimension.first != 0 && n != m_dimension.first) {
        throw std::invalid_argument("x-dimension is invalid");
    }
    double x0_2 = x[0] / 2.0;
    for (size_t k = 0; k < n; k++) {
        double yk = x0_2;
        for (size_t i = 1; i < n; i++) {
            yk += (x[i] * std::cos(i * M_PI * (k + 0.5) / n));
        }
        y.set(k, 0, gamma * y[k] + alpha * yk);
    }
    return ForBESUtils::STATUS_OK;
}

int OpDCT3::callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma) {
    size_t n = x.length();
    if (m_dimension.first != 0 && n != m_dimension.first) {
        throw std::invalid_argument("x-dimension is invalid");
    }
    double tk = 0.0;
    for (size_t i = 0; i < n; i++) {
        tk += x[i] / 2.0;
    }
    y.set(0, 0, gamma * y.get(0, 0) + alpha * tk);
    for (size_t k = 1; k < n; k++) {
        tk = 0.0;
        for (size_t i = 0; i < n; i++) {
            tk += (x[i] * std::cos(k * M_PI * (i + 0.5) / n));
        }
        y.set(k, 0, gamma * y[k] + alpha * tk);
    }
    return ForBESUtils::STATUS_OK;
}

std::pair<size_t, size_t> OpDCT3::dimensionIn() {
    return m_dimension;
}

std::pair<size_t, size_t> OpDCT3::dimensionOut() {
    return m_dimension;
}

bool OpDCT3::isSelfAdjoint() {
    return false;
}






