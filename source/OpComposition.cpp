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

OpComposition::~OpComposition() {
}

Matrix OpComposition::call(Matrix& x) {
    Matrix y = B.call(x);   // y = B(x)
    return A.call(y);       // z = A(y)
}

size_t OpComposition::dimensionIn() {
    return B.dimensionIn();
}

size_t OpComposition::dimensionOut() {
    return A.dimensionOut();
}

bool OpComposition::isSelfAdjoint() {
    return A.isSelfAdjoint() && B.isSelfAdjoint();
}

Matrix OpComposition::callAdjoint(Matrix& x) {
    Matrix y = A.callAdjoint(x); // y = B(x)
    return B.callAdjoint(y); // z = A(y)
}
