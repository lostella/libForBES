/* 
 * File:   OpDCT2.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 15, 2015, 3:36 PM
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

#include "OpDCT2.h"

OpDCT2::OpDCT2() : LinearOperator(), m_dimension(_EMPTY_OP_DIM) {

}

OpDCT2::OpDCT2(size_t n) : m_dimension(_VECTOR_OP_DIM(n)) {
}

OpDCT2::~OpDCT2() {
}

int OpDCT2::call(Matrix& y, double alpha, Matrix& x, double gamma) {
    size_t n = x.getNrows();
    for (size_t k = 0; k < n; k++) {
        double yk = 0.0;
        double aik;
        for (size_t i = 0; i < n; i++) {
            aik = std::cos(M_PI * (static_cast<double> (i) + 0.5) * static_cast<double> (k) / static_cast<double> (n));
            yk += x.get(i, 0) * aik;
        }
        y.set(k, 0, gamma * y.get(k, 0) + alpha * yk);
    }
    return ForBESUtils::STATUS_OK;
}

int OpDCT2::callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma) {
    size_t n = x.getNrows();
    for (size_t k = 0; k < n; k++) {
        double v = 0.0;
        double aki;
        for (size_t i = 0; i < n; i++) {
            aki = std::cos(M_PI * (static_cast<double> (k) + 0.5) * static_cast<double> (i) / static_cast<double> (n));
            v += x.get(i, 0) * aki;
        }
        y.set(k, 0, gamma * y.get(k, 0) + alpha * v);
    }
    return ForBESUtils::STATUS_OK;
}

std::pair<size_t, size_t> OpDCT2::dimensionIn() {
    return m_dimension;
}

std::pair<size_t, size_t> OpDCT2::dimensionOut() {
    return m_dimension;
}

bool OpDCT2::isSelfAdjoint() {
    return false;
}
