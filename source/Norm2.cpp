/* 
 * File:   Norm2.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on October 30, 2015, 6:54 PM
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

#include "Norm2.h"
#include <cmath>

double vecNorm2(Matrix& x);

Norm2::Norm2() : Norm() {
    m_mu = 1.0;
}

Norm2::Norm2(double mu) : Norm(), m_mu(mu) {
}

Norm2::~Norm2() {
}

inline double vecNorm2(Matrix& x) {
    double norm_x = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        norm_x += std::pow(x[i], 2);
    }
    norm_x = std::sqrt(norm_x);
    return norm_x;
}

int Norm2::call(Matrix& x, double& f) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {        
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    f = m_mu * vecNorm2(x);
    return ForBESUtils::STATUS_OK;
}

int Norm2::dualNorm(Matrix& x, double& norm) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }   
    //LCOV_EXCL_STOP
    norm = vecNorm2(x) / m_mu;
    return ForBESUtils::STATUS_OK;
}

int Norm2::callProx(Matrix& x, double gamma, Matrix& prox) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    double norm_x = vecNorm2(x);
    double gm = gamma*m_mu;
    if (norm_x > gm) {
        double s = 1 - gm / norm_x;
        prox = s * x;
    }
    return ForBESUtils::STATUS_OK;
}

int Norm2::callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    double norm_x = vecNorm2(x);    
    double gm = gamma*m_mu;
    if (norm_x > gm) {
        double s = 1 - gm / norm_x;
        for (size_t i = 0; i < x.getNrows(); i++) {
            prox.set(i, 0, s * x[i]); // prox = s * x;
        }
        f_at_prox = m_mu * s * norm_x;
    } else {
        f_at_prox = 0.0;
    }
    return ForBESUtils::STATUS_OK;
}

FunctionOntologicalClass Norm2::category() {
    return FunctionOntologyRegistry::norm();
}






