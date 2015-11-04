/* 
 * File:   OpLinearCombination.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 14, 2015, 9:25 PM
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

#include "OpLinearCombination.h"

OpLinearCombination::OpLinearCombination(LinearOperator& A, LinearOperator& B, double a, double b) :
LinearOperator(), m_A(A), m_B(B), m_a(a), m_b(b) {
    if (A.dimensionIn() != B.dimensionIn()) {
        throw std::invalid_argument("A and B have incompatible input dimensions");
    }
    if (A.dimensionOut() != B.dimensionOut()) {
        throw std::invalid_argument("A and B have incompatible input dimensions");
    }
}

OpLinearCombination::~OpLinearCombination() {
}

