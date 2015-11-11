/* 
 * File:   OpComposition.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 14, 2015, 9:26 PM
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

#include "OpComposition.h"

OpComposition::OpComposition(LinearOperator& A, LinearOperator& B) : LinearOperator(), m_A(A), m_B(B) {
    // check dimensions
    if (A.dimensionIn() != B.dimensionOut()) {
        throw std::invalid_argument("A and B have incompatible dimensions; AoB is not well defined.");
    }
}

OpComposition::~OpComposition() {
}

int OpComposition::call(Matrix& y, double alpha, Matrix& x, double gamma) {
    Matrix t(m_B.dimensionOut().first, m_B.dimensionOut().second);
    int status = m_B.call(t, 1.0, x, 0.0); // t = B(x)
    if (ForBESUtils::is_status_error(status)) {
        return status;
    }
    status = std::max(status, m_A.call(y, alpha, t, gamma));
    return status;
}

int OpComposition::callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma) {
    Matrix t = m_A.callAdjoint(x);            
    return m_B.callAdjoint(y, alpha, t, gamma);
}

std::pair<size_t, size_t> OpComposition::dimensionIn() {
    return m_B.dimensionIn();
}

std::pair<size_t, size_t> OpComposition::dimensionOut() {
    return m_A.dimensionOut();
}

bool OpComposition::isSelfAdjoint() {
    return m_A.isSelfAdjoint() && m_B.isSelfAdjoint();
}


