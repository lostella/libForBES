/* 
 * File:   QuadraticOperator.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on July 24, 2015, 8:49 PM
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

#include "QuadraticOperator.h"
#include "ForBESUtils.h"

QuadraticOperator::QuadraticOperator(LinearOperator& T) :
Function(), T(T) {
    if (T.dimensionIn() != T.dimensionOut()) {
        throw std::invalid_argument("T has incompatible dimensions");
    }
}

QuadraticOperator::~QuadraticOperator() {
}

int QuadraticOperator::call(Matrix& x, double& f, Matrix& grad) {
    int status = computeGradient(x, grad);
    f = (x * grad).get(0, 0)/2;
    return status;
}

int QuadraticOperator::call(Matrix& x, double& f) {
    Matrix grad;
    int status = call(x, f, grad);
    return status;
}

int QuadraticOperator::computeGradient(Matrix& x, Matrix& grad) {
    grad = T.call(x); /* Tx = T[x]; assertion: Tx is a vector of the same length as x */
    
    if (grad.isEmpty()) {
        throw std::logic_error("Linear operator returned an empty vector");
    }
    if (!grad.isColumnVector()) {
        throw std::logic_error("Linear operator did not return a column vector");
    }
    if (grad.getNrows() != x.getNrows()) {
        throw std::logic_error("Incompatible dimensions between T(x) and x");
    }
    return ForBESUtils::STATUS_OK;
}

