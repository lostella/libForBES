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
#include <sstream>

void call_1d(Matrix & Tx, Matrix& x, const size_t n, double alpha);
void call_1d(Matrix & Tx, Matrix& x, const size_t n, double alpha, double gamma);
void callAdjoint_1d(Matrix& Tstar_x, Matrix& y, const size_t n, double alpha);
void callAdjoint_1d(Matrix& Tstar_x, Matrix& y, const size_t n, double alpha, double gammna);

OpGradient::OpGradient() : m_dimension(_EMPTY_OP_DIM) {
}

OpGradient::OpGradient(size_t n) : LinearOperator(), m_dimension(_VECTOR_OP_DIM(n)) {

}

OpGradient::~OpGradient() {
}

/**
 * Calls OpGradient when the input is 1D
 * @param Tx matrix (vector) to store the result
 * @param x input vector x
 * @param n size of x
 */
void call_1d(Matrix & Tx, Matrix& x, const size_t n, double alpha) {
    for (size_t i = 0; i < n - 1; i++) {
        Tx[i] = alpha * (x[i + 1] - x[i]);
    }
}

void call_1d(Matrix & Tx, Matrix& x, const size_t n, double alpha, double gamma) {
    for (size_t i = 0; i < n - 1; i++) {
        Tx[i] = gamma * Tx[i] + alpha * (x[i + 1] - x[i]);
    }
}

void callAdjoint_1d(Matrix& Tstar_x, Matrix& y, const size_t n, double alpha) {
    Tstar_x[0] = -alpha * y[0];
    for (size_t i = 1; i < n - 1; i++) {
        Tstar_x[i] = alpha * (y[i - 1] - y[i]);
    }
    Tstar_x[n - 1] = alpha * y[n - 2];
}

void callAdjoint_1d(Matrix& Tstar_x, Matrix& y, const size_t n, double alpha, double gamma) {
    Tstar_x.set(0, 0, gamma * Tstar_x.get(0, 0) - alpha * y.get(0, 0));
    for (size_t i = 1; i < n - 1; i++) {
        Tstar_x[i] = gamma * Tstar_x[i] + alpha * (y[i - 1] - y[i]);
    }
    Tstar_x[n - 1] = gamma * Tstar_x[n - 1] + alpha * y[n - 2];
}

int OpGradient::call(Matrix& y, double alpha, Matrix& x, double gamma) {
    const size_t n = x.getNrows();
    if (m_dimension.first != 0 && n != m_dimension.first) {
        std::ostringstream oss;
        oss << "[call] OpGradient operator with dimension " << m_dimension.first
                << "; argument is of incompatible dimensions " << n << "x" << x.getNcols();
        throw std::invalid_argument(oss.str().c_str());
    }
    if (n <= 1) {
        y = alpha*x;
        return ForBESUtils::STATUS_NUMERICAL_PROBLEMS;
    }
    call_1d(y, x, n, alpha);
    return ForBESUtils::STATUS_OK;
}

int OpGradient::callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma) {
    callAdjoint_1d(y, x, x.getNrows() + 1, alpha);
    return ForBESUtils::STATUS_OK;
}

std::pair<size_t, size_t> OpGradient::dimensionIn() {
    return m_dimension;
}

std::pair<size_t, size_t> OpGradient::dimensionOut() {
    std::pair<size_t, size_t> dims;
    dims.second = static_cast<size_t> (1);
    if (m_dimension.first == 0) {
        dims.first = static_cast<size_t> (0);
    } else {
        dims.first = static_cast<size_t> (dimensionIn().first - 1);
    }
    return dims;
}

bool OpGradient::isSelfAdjoint() {
    return false;
}






