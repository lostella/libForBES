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
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    f = 0.0;
    double xi = 0.0;
    for (size_t i = 0; i < x.getNrows(); i++) {
        xi = x.get(i, 0);
        f += mu * std::abs(xi) + (lambda / 2.0) * std::pow(xi, 2);
    }
    return ForBESUtils::STATUS_OK;
}

int ElasticNet::callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    double xi = 0.0;
    double gm = gamma*mu;
    double alpha = (1+lambda*gamma);
    double tau;
    for (size_t i = 0; i < x.getNrows(); i++) {
        xi = x.get(i, 0);
        if  (xi < 0) {
            prox.set()
            tau = -max(0.0, abs(xi) - gm)/alpha;
        } else {
            tau = max(0.0, abs(xi) - gm)/alpha;
        }
    }
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

int ElasticNet::callProx(const Matrix& x, double gamma, Matrix& prox) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

int ElasticNet::category() {
    return 0;
}