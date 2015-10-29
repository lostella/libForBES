/* 
 * File:   ElasticNet.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on October 28, 2015, 7:43 PM
 */

#include "ElasticNet.h"

ElasticNet::ElasticNet(double lambda, double mu) : Function(), lambda(lambda), mu(mu) {
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
    double xi = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        xi = x.get(i, 0);
        f += mu * std::abs(xi) + (lambda / 2.0) * std::pow(xi, 2);
    }
    return ForBESUtils::STATUS_OK;
}

int ElasticNet::callProx(const Matrix& x, double gamma, Matrix& prox, double& g_at_prox) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    double xi = 0.0;
    double gm = gamma * mu;
    double alpha = 1 + lambda * gamma; // alpha > 0 [assuming gamma>0 and lambda>0].
    double yi = 0.0;
    g_at_prox = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        xi = x.get(i, 0);
        yi = max(0.0, abs(xi) - gm) / alpha;
        prox.set(i, 0, (xi < 0 ? -1 : 1) * yi);
        g_at_prox += mu * yi + (lambda / 2.0) * std::pow(yi, 2);
    }
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

int ElasticNet::callProx(const Matrix& x, double gamma, Matrix& prox) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP
    double xi = 0.0;
    double gm = gamma * mu;
    double alpha = 1 + lambda * gamma; // alpha > 0 [assuming gamma>0 and lambda>0].
    for (size_t i = 0; i < x.getNrows(); i++) {
        xi = x.get(i, 0);
        prox.set(i, 0, (xi < 0 ? -1 : 1) * max(0.0, abs(xi) - gm) / alpha);
    }
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}
