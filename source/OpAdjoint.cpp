/* 
 * File:   OpAdjoint.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 15, 2015, 2:57 PM
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

#include "OpAdjoint.h"

OpAdjoint::OpAdjoint(LinearOperator& op) : LinearOperator(), m_originalOperator(op) {

}

OpAdjoint::~OpAdjoint() {
}

Matrix OpAdjoint::call(Matrix& x) {
    return m_originalOperator.callAdjoint(x);
}

Matrix OpAdjoint::callAdjoint(Matrix& x) {
    return m_originalOperator.call(x);
}

std::pair<size_t, size_t> OpAdjoint::dimensionIn() {
    return m_originalOperator.dimensionOut();
}

std::pair<size_t, size_t> OpAdjoint::dimensionOut() {    
    return m_originalOperator.dimensionIn();
}

bool OpAdjoint::isSelfAdjoint() {
    return m_originalOperator.isSelfAdjoint();
}






