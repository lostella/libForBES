/* 
 * File:   QuadOverAffine.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on July 24, 2015, 4:55 PM
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

#include "QuadOverAffine.h"


QuadOverAffine::QuadOverAffine() {
}

QuadOverAffine::QuadOverAffine(const QuadOverAffine& orig) {
}

QuadOverAffine::~QuadOverAffine() {
}

QuadOverAffine::QuadOverAffine(Matrix& Q, Matrix& q, Matrix& A, Matrix& b) {
    this->Q = &Q;
    this->q = &q;
    this->A = &A;
    this->b = &b;    
}


