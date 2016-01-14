/* 
 * File:   OpLTI.cpp
 * Author: chung
 * 
 * Created on September 30, 2015, 6:20 PM
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

#include "OpLTI.h"
#include "MatrixFactory.h"
#include <sstream>
#include <iostream>


OpLTI::OpLTI(Matrix& A, Matrix& B, size_t N_horizon) :
LinearOperator(),
m_A(A),
m_B(B),
m_N(N_horizon) {
    if (A.getNrows() != A.getNcols()) {
        throw std::invalid_argument("System matrix A must be square");
    }
    if (A.getNrows() != B.getNrows()) {
        throw std::invalid_argument("A and B have incompatible dimensions");
    }
    if (N_horizon == 0) {
        throw std::invalid_argument("N_horizon cannot be zero");
    }
}

OpLTI::~OpLTI() {
}

int OpLTI::call(Matrix& y, double alpha, Matrix& u, double gamma) {
    // y := gamma * y + alpha * T(u)
    // T(u) = [x1, x2, ..., xN]
    if (!u.isColumnVector()) {
        throw std::invalid_argument("u should be a column vector");
    }
    if (u.getNrows() != (m_N - 1) * m_B.getNcols()) {
        std::ostringstream oss;
        oss << "u: wrong dimension - should be (N-1)*n_u = "
                << (m_N - 1) << "*" << m_B.getNcols()
                << "=" << ((m_N - 1) * m_B.getNcols());
        throw std::invalid_argument(oss.str().c_str());
    }
    size_t n = m_A.getNrows();
    size_t m = m_B.getNcols();
    const size_t zero = 0;
    Matrix y0 = MatrixFactory::ShallowVector(y, n, zero);
    Matrix u0 = MatrixFactory::ShallowVector(u, m, zero);
    int status = Matrix::mult(y0, alpha, m_B, u0, gamma);
    ForBESUtils::fail_on_error(status);
    std::cout << std::endl;
    for (size_t i = 1; i < m_N - 1; i++) {
        Matrix ui = MatrixFactory::ShallowVector(u, m, i * m);
        std::cout << "u[" << i << "] = " << ui << std::endl;
    }
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

int OpLTI::callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

std::pair<size_t, size_t> OpLTI::dimensionIn() {
    std::pair<size_t, size_t> dims(static_cast<size_t> (0), static_cast<size_t> (1));
    return dims;
}

std::pair<size_t, size_t> OpLTI::dimensionOut() {
    std::pair<size_t, size_t> dims(static_cast<size_t> (0), static_cast<size_t> (1));
    return dims;
}

bool OpLTI::isSelfAdjoint() {
    return false;
}


