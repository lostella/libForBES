/* 
 * File:   ForBESUtils.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on July 24, 2015, 5:17 PM
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

#include <stdexcept>

#include "ForBESUtils.h"

/*
 * Status codes 0~100    are informative
 *              200~400  are warnings
 *              500~1000 are errors
 */

const int ForBESUtils::STATUS_OK = 0;
const int ForBESUtils::STATUS_HAD_TO_REALLOC = 1;

const int ForBESUtils::STATUS_NUMERICAL_PROBLEMS = 500;
const int ForBESUtils::STATUS_UNDEFINED_FUNCTION = 501;

void ForBESUtils::fail_on_error(int status) {
    if (is_status_error(status)) {
        throw std::logic_error("Failed");
    }
}

bool ForBESUtils::is_status_error(int status) {
    return (status > _FORBES_ERROR_MIN && status < _FORBES_ERROR_MAX);
}

bool ForBESUtils::is_status_ok(int status) {
    return (status > _FORBES_OK_MIN && status < _FORBES_OK_MAX);
}

bool ForBESUtils::is_status_warning(int status) {
    return (status > _FORBES_WARNING_MIN && status < _FORBES_WARNING_MAX);
}








