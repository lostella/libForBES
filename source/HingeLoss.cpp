/* 
 * File:   HingeLoss.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on October 29, 2015, 10:49 PM
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

#include "HingeLoss.h"

HingeLoss::HingeLoss() : Function() {
    m_mu = 1.0;
}

HingeLoss::HingeLoss(Matrix* b, double mu) :
Function(), m_b(b), m_mu(mu) {
}

HingeLoss::HingeLoss(Matrix* b) :
Function(), m_b(b) {
    m_mu = 1.0;
}

HingeLoss::~HingeLoss() {
}

int HingeLoss::call(Matrix& x, double& f) {
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    f = 0.0;
    double si = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        si = 1 - m_b->get(i, 0) * x.get(i, 0);
        if (si > 0) {
            f += si;
        }
    }
    f /= m_mu;
    return ForBESUtils::STATUS_OK;
}

int HingeLoss::callProx(const Matrix& x, double gamma, Matrix& prox) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

int HingeLoss::callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

FunctionOntologicalClass HingeLoss::category() {
    return FunctionOntologyRegistry::loss();
}
