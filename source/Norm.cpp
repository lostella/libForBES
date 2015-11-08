/* 
 * File:   Norm.cpp
 * Author: chung
 * 
 * Created on October 30, 2015, 5:32 PM
 */

#include "Norm.h"

Norm::Norm() : Function() {
}

Norm::~Norm() {
}

int Norm::callConj(Matrix& x, double& f_star) {
    double dNorm;
    int status = dualNorm(x, dNorm);
    if (ForBESUtils::STATUS_OK != status) {
        return status;
    }
    if (dNorm <= 1) {
        f_star = 0.0;
    } else {
        f_star = INFINITY;
    }
    return ForBESUtils::STATUS_OK;
}

//LCOV_EXCL_START
int Norm::dualNorm(Matrix& x, double& norm) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

FunctionOntologicalClass Norm::category() {
    return FunctionOntologyRegistry::norm();
}
//LCOV_EXCL_STOP