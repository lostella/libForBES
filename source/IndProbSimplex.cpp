/* 
 * File:   IndProbSimplex.cpp
 * Author: chung
 * 
 * Created on January 15, 2016, 2:37 AM
 */

#include "IndProbSimplex.h"
#include <vector>
#include <algorithm>
#include <cmath>

IndProbSimplex::IndProbSimplex() {
}

IndProbSimplex::~IndProbSimplex() {
}

/**
 * Computes the value of g(t;x) = 1'(x-t1)-1
 * @param x reference to input vector
 * @param t current estimate of t
 * @return g(t;x)
 */
double compute_f893301(Matrix &x, double t) {
    double g = 0.0;
    double diff;
    for (size_t i = 0; i < x.getNrows(); i++) {
        diff = x[i] - t;
        if (diff > 0) {
            g += diff;
        }
    }
    return g - 1.0;
}

int IndProbSimplex::callProx(Matrix& x, double gamma, Matrix& prox) {
    std::vector<double> Avec(x.getData(), x.getData() + x.getNrows());
    double max_val = *std::max_element(Avec.begin(), Avec.end());

    /* BISECTION */
    double t0 = max_val - 1.0;
    double t1 = max_val;
    double t2;
    double f2 = 1;
    const size_t max_iter = 1000;
    const double tol = 1e-10;
    size_t i = 0;
    do {
        t2 = (t0 + t1) / 2;
        double f0 = compute_f893301(x, t0);
        double f1 = compute_f893301(x, t1);
        f2 = compute_f893301(x, t2);
        if (f0 * f2 < 0) {
            t1 = t2;
        } else {
            t0 = t2;
        }
        i++;
    } while (std::abs(f2) > tol && i < max_iter);

    /* CHECK QUALITY OF SOLUTION (max_iter, tolerance) */
    if (i == max_iter && std::abs(f2) > tol) {
        return ForBESUtils::STATUS_NUMERICAL_PROBLEMS;
    }

    /* CONSTRUCT THE PROX */
    for (size_t i = 0; i < prox.getNrows(); i++) {
        prox[i] = std::max(0.0, x[i] - t2);
    }

    return ForBESUtils::STATUS_OK;
}

int IndProbSimplex::callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    f_at_prox = 0.0;
    return callProx(x, gamma, prox);
}

int IndProbSimplex::call(Matrix& x, double& f) {
    /* Implement me! */
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

FunctionOntologicalClass IndProbSimplex::category() {
    return FunctionOntologyRegistry::indicator();
}
