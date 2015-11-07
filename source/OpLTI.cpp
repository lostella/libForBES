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

OpLTI::OpLTI(Matrix& A, Matrix& B) : LinearOperator(), A(A), B(B) {
    if (A.getNrows() != B.getNrows()){
        throw std::invalid_argument("A and B have incompatible dimensions");
    }
}

OpLTI::~OpLTI() {
}

Matrix OpLTI::call(Matrix& u) {
    size_t n = A.getNrows();
    size_t m = B.getNcols();
    
    if (u.isColumnVector()){
        
    } else {
        
    }    
    throw std::logic_error("NIY");
}

Matrix OpLTI::callAdjoint(Matrix& x) {
    throw std::logic_error("NIY");
}

std::pair<size_t, size_t> OpLTI::dimensionIn() {
    std::pair<size_t, size_t> dims(static_cast<size_t>(0),static_cast<size_t>(1));
    return dims;
}

std::pair<size_t, size_t> OpLTI::dimensionOut() {
    std::pair<size_t, size_t> dims(static_cast<size_t>(0),static_cast<size_t>(1));
    return dims;
}

bool OpLTI::isSelfAdjoint() {
    return false;
}


