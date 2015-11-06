/* 
 * File:   Function.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on July 9, 2015, 3:35 AM
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

#include "Function.h"

Function::Function() { }

Function::~Function() {
}

//LCOV_EXCL_START
int Function::call(Matrix& x, double& f, Matrix& grad)  {
    int status;
    status = computeGradient(x, grad);
    if (ForBESUtils::STATUS_OK != status) {
        return status;
    }
    status = call(x, f);
    if (ForBESUtils::STATUS_OK != status) {
        return status;
    }    
    return ForBESUtils::STATUS_OK;
}

int Function::call(Matrix& x, double& f) {    
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}
//LCOV_EXCL_STOP

int Function::callConj(const Matrix& x, double& f_star) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

int Function::callConj(const Matrix& x, double& f_star, Matrix& grad) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

int Function::callProx(const Matrix& x, double gamma, Matrix& prox) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

//LCOV_EXCL_START
int Function::callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}


int Function::computeGradient(Matrix& x, Matrix& grad) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}
//LCOV_EXCL_STOP
