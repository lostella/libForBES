/* 
 * File:   LinearOperator.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on July 24, 2015, 8:44 PM
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

#include "LinearOperator.h"

LinearOperator::LinearOperator() {
}

LinearOperator::~LinearOperator() {
}

Matrix LinearOperator::call(Matrix& x) {
    Matrix y(dimensionOut().first, dimensionOut().second);
    const double gamma = 0.0;
    const double alpha = 1.0;
    ForBESUtils::fail_on_error(call(y, alpha, x, gamma));
    return y;
}

Matrix LinearOperator::callAdjoint(Matrix& x) {
    Matrix y_star(dimensionIn().first, dimensionIn().second);
    const double gamma = 0.0;
    const double alpha = 1.0;
    ForBESUtils::fail_on_error(callAdjoint(y_star, alpha, x, gamma));
    return y_star;
}
