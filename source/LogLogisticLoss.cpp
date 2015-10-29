/* 
 * File:   LogLogisticLoss.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on October 29, 2015, 5:08 PM
 */

#include "LogLogisticLoss.h"

LogLogisticLoss::LogLogisticLoss() {
    mu = 1.0;
}

LogLogisticLoss::~LogLogisticLoss() {
}

LogLogisticLoss::LogLogisticLoss(double mu) :
Function(), mu(mu) {
}

int LogLogisticLoss::call(Matrix& x, double& f, Matrix& grad) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    f = 0.0;
    long double si = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        si = std::exp(x.get(i, 0));
        si = si / (1 + si);
        f -= std::log(si);
        grad.set(i, 0, mu * (si - 1));
    }
    f *= mu;
    return ForBESUtils::STATUS_OK;
}

int LogLogisticLoss::call(Matrix& x, double& f) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    f = 0.0;
    double si = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        si = std::exp(x.get(i, 0));
        si = si / (1 + si);
        f -= std::log(si);
    }
    f *= mu;
    return ForBESUtils::STATUS_OK;
}

FunctionOntologicalClass LogLogisticLoss::category() {
    return FunctionOntologyRegistry::loss();
}





