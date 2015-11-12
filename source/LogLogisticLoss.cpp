/* 
 * File:   LogLogisticLoss.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on October 29, 2015, 5:08 PM
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

#include "LogLogisticLoss.h"
#include <cmath>

LogLogisticLoss::LogLogisticLoss() {
    m_mu = 1.0;
}

LogLogisticLoss::~LogLogisticLoss() {
}

LogLogisticLoss::LogLogisticLoss(double mu) :
Function(), m_mu(mu) {
}

int LogLogisticLoss::call(Matrix& x, double& f, Matrix& grad) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    f = 0.0;
    int status = ForBESUtils::STATUS_OK;
    for (size_t i = 0; i < x.getNrows(); i++) {        
        double xi = x.get(i,0);
        if (xi < 33) {                      /* because for values higher than 33, 
                                             * s1 is practically equal to 1. 
                                             * This saves the computational burden for
                                             * high values of x_i */
            double si = std::exp(xi);       /* si = e^xi              */
            si = si / (1 + si);             /* si = e^xi / (1+e^xi)   */
            f -= std::log(si);              /* f -= ln(si)            */
            grad.set(i, 0, m_mu * (si - 1));  /* prox_i = mu*(si-1)     */
        }
    }
    f *= m_mu;
    return status;
}

int LogLogisticLoss::call(Matrix& x, double& f) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    f = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        double si = std::exp(x.get(i, 0));
        si /= (1.0 + si);
        f -= std::log(si);
    }
    f *= m_mu;
    return ForBESUtils::STATUS_OK;
}

FunctionOntologicalClass LogLogisticLoss::category() {
    FunctionOntologicalClass logLogisticLoss("LogLogisticLoss");
    logLogisticLoss.set_defines_f(true);
    logLogisticLoss.set_defines_grad(true);
    logLogisticLoss.getSuperclasses().push_back(FunctionOntologyRegistry::loss());
    return logLogisticLoss;
}





