/* 
 * File:   OpGradient.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 16, 2015, 1:35 AM
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

#include "OpGradient.h"

OpGradient::OpGradient() {
    m_dimension = 0;
}

OpGradient::OpGradient(size_t n) : LinearOperator(), m_dimension(n) {
}

OpGradient::~OpGradient() {
}

/**
 * Calls OpGradient when the input is 1D
 * @param Tx matrix (vector) to store the result
 * @param x input vector x
 * @param n size of x
 */
void call_1d(Matrix & Tx, Matrix& x, const size_t n) {
    for (size_t i = 0; i < n - 1; i++) {
        Tx.set(i, 0, x.get(i + 1, 0) - x.get(i, 0));
    }
}

void callAdjoint_1d(Matrix& Tstar_x, Matrix& y, const size_t n) {
    Tstar_x.set(0, 0, -y.get(0, 0));
    for (size_t i = 1; i < n - 1; i++) {
        Tstar_x.set(i, 0, y.get(i - 1, 0) - y.get(i, 0));
    }
    Tstar_x.set(n - 1, 0, y.get(n - 2, 0));
}

Matrix OpGradient::call(Matrix& x) {
    const size_t n = x.length();
    if (m_dimension != 0 && n != m_dimension) {
        throw std::invalid_argument("x-dimension is invalid");
    }
    Matrix Tx(n - 1, 1);
    if (n <= 1) {
        return Tx;
    }
    call_1d(Tx, x, n);
    return Tx;
}

Matrix OpGradient::callAdjoint(Matrix& y) {
    const size_t n = y.length() + 1;
    if (m_dimension != 0 && n != m_dimension) {
        throw std::invalid_argument("x-dimension is invalid");
    }
    Matrix Tstar_x(n, 1);
    callAdjoint_1d(Tstar_x, y, n);
    return Tstar_x;
}

size_t OpGradient::dimensionIn() {
    return m_dimension;
}

size_t OpGradient::dimensionOut() {
    return m_dimension == 0 ? 0 : (m_dimension - 1);
}

bool OpGradient::isSelfAdjoint() {
    return false;
}






