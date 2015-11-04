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

Matrix OpComposition::call(Matrix& x) {
    Matrix y = m_B.call(x);   // y = B(x)
    return m_A.call(y);       // z = A(y)
}

size_t OpComposition::dimensionIn() {
    return m_B.dimensionIn();
}

size_t OpComposition::dimensionOut() {
    return m_A.dimensionOut();
}

bool OpComposition::isSelfAdjoint() {
    return m_A.isSelfAdjoint() && m_B.isSelfAdjoint();
}

Matrix OpComposition::callAdjoint(Matrix& x) {
    /*
     * T(x) = A(B(x))
     * so
     * T*(x) = B*(A*(x)) = B*(y), where y = A*(x)
     */
    Matrix y = m_A.callAdjoint(x); // y = A*(x)
    return m_B.callAdjoint(y);     // z = B*(y)
}
