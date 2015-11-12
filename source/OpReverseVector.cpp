/* 
 * File:   OpReverseVector.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 15, 2015, 12:57 PM
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

#include "OpReverseVector.h"
#include <algorithm>
#include <iterator>
#include <limits>
#include <cmath>

OpReverseVector::OpReverseVector() {
    m_vectorDim = 0;
}

OpReverseVector::~OpReverseVector() {
}

OpReverseVector::OpReverseVector(size_t n) {
    this->m_vectorDim = n;
}

int OpReverseVector::call(Matrix& y, double alpha, Matrix& x, double gamma) {
    if (y.getNrows() == 0) {
        y = Matrix(x.getNrows(), 1);
    }
    bool is_gamma_zero = (std::abs(gamma) < std::numeric_limits<double>::epsilon());
    size_t n = x.length();
    for (size_t i = 0; i < n / 2; ++i) {
        double temp;
        temp = x.get(n - i - 1, 0);
        if (is_gamma_zero) {
            y.set(n - i - 1, 0, alpha * x.get(i, 0));
            y.set(i, 0, alpha * temp);
        } else {
            y.set(n - i - 1, 0, gamma * y.get(n - i - 1, 0) + alpha * x.get(i, 0));
            y.set(i, 0, gamma * y.get(i, 0) + alpha * temp);
        }
    }
    if (n % 2 == 1) {
        size_t middle_idx = n / 2;
        if (is_gamma_zero) {
            y.set(middle_idx, 0, alpha * x.get(middle_idx, 0));
        } else {
            y.set(middle_idx, 0, gamma * y.get(middle_idx, 0) + alpha * x.get(middle_idx, 0));
        }
    }
    return ForBESUtils::STATUS_OK;
}

int OpReverseVector::callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma) {
    return call(y, alpha, x, gamma);
}

std::pair<size_t, size_t> OpReverseVector::dimensionIn() {
    return _VECTOR_OP_DIM(m_vectorDim);
}

std::pair<size_t, size_t> OpReverseVector::dimensionOut() {
    return _VECTOR_OP_DIM(m_vectorDim);
}

bool OpReverseVector::isSelfAdjoint() {
    return true;
}







