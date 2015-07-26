/* 
 * File:   Function.cpp
 * Author: chung
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

const int Function::CAT_QUADRATIC = 100;

const int Function::CAT_INDICATOR = 150;

const int Function::CAT_UNCATEGORIZED = 9999;

Function::Function() {

}

Function::Function(const Function& orig) {
}

Function::~Function() {
}

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
