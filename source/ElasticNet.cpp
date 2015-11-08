/* 
 * File:   ElasticNet.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on October 28, 2015, 7:43 PM
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

#include "ElasticNet.h"

ElasticNet::ElasticNet(double lambda, double mu) : Function(), m_mu(mu), m_lambda(lambda) {
}

ElasticNet::~ElasticNet() {
}

int ElasticNet::call(Matrix& x, double& f) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    f = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        double xi;
        xi = x.get(i, 0);
        f += m_mu * std::abs(xi) + (m_lambda / 2.0) * std::pow(xi, 2);
    }
    return ForBESUtils::STATUS_OK;
}

int ElasticNet::callProx(Matrix& x, double gamma, Matrix& prox, double& g_at_prox) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    double gm = gamma * m_mu;
    double alpha = 1 + m_lambda * gamma; // alpha > 0 [assuming gamma>0 and lambda>0].
    g_at_prox = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        double xi;
        double yi;
        xi = x.get(i, 0);
        yi = max(0.0, abs(xi) - gm) / alpha;
        prox.set(i, 0, (xi < 0 ? -1 : 1) * yi);
        g_at_prox += m_mu * yi + (m_lambda / 2.0) * std::pow(yi, 2);
    }
    return ForBESUtils::STATUS_OK;
}

int ElasticNet::callProx(Matrix& x, double gamma, Matrix& prox) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    double gm = gamma * m_mu;
    double alpha = 1 + m_lambda * gamma; // alpha > 0 [assuming gamma>0 and lambda>0].
    for (size_t i = 0; i < x.getNrows(); i++) {
        double xi;
        xi = x.get(i, 0);
        prox.set(i, 0, (xi < 0 ? -1 : 1) * max(0.0, abs(xi) - gm) / alpha);
    }
    return ForBESUtils::STATUS_OK;
}

FunctionOntologicalClass ElasticNet::category() {
    FunctionOntologicalClass ont("ElasticNet");
    ont.set_defines_conjugate(false);
    ont.set_defines_conjugate_grad(false);
    ont.set_defines_f(true);
    ont.set_defines_grad(false);
    ont.set_defines_prox(true);
    ont.getSuperclasses().push_back(FunctionOntologyRegistry::function());
    return ont;
}
