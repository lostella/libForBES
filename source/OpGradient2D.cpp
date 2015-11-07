/* 
 * File:   OpGradient2D.cpp
 * Author: chung
 * 
 * Created on September 16, 2015, 6:20 PM
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

#include "OpGradient2D.h"

OpGradient2D::OpGradient2D() {
}

OpGradient2D::~OpGradient2D() {
}

Matrix OpGradient2D::call(Matrix& x) {
    throw std::logic_error("NIY");
}

Matrix OpGradient2D::callAdjoint(Matrix& x) {
    throw std::logic_error("NIY");
}

std::pair<size_t, size_t> OpGradient2D::dimensionIn() {
    throw std::logic_error("NIY");
}

std::pair<size_t, size_t> OpGradient2D::dimensionOut() {
    throw std::logic_error("NIY");
}

bool OpGradient2D::isSelfAdjoint() {
    throw std::logic_error("NIY");
}

