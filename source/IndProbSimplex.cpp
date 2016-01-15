/* 
 * File:   IndProbSimplex.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on January 15, 2016, 2:37 AM
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

#include "IndProbSimplex.h"
#include <vector>
#include <algorithm>
#include <cmath>

IndProbSimplex::IndProbSimplex() {

}

IndProbSimplex::~IndProbSimplex() {
}

int IndProbSimplex::callProx(Matrix& x, double gamma, Matrix& prox) {

    size_t n = x.getNrows();

    /* x_hat := x */
    Matrix x_hat(x);
    for (size_t j = 0; j < n; j++) {
        if (x_hat[j] < 0.0) {
            x_hat[j] = 0.0;
        }
    }

    /* x_hat := rev_sort(x_hat) */
    std::vector<double> x_hat_vec(x_hat.getData(), x_hat.getData() + n);
    //std::sort(x_hat_vec.begin(), x_hat_vec.end(), std::greater<double>());
    std::sort(x_hat_vec.rbegin(), x_hat_vec.rend());

    double t = 0.0;
    bool flag = true;
    size_t i = 0;
    double val;
    do {
        val = x_hat_vec[i];
        t += val;
        flag = val * (i + 1.0) > t - 1.0;
        ++i;
    } while (flag && i < n - 1);
    t -= 1.0 + val;

    double theta = std::max(0.0, t / (i - 1));

    for (size_t j = 0; j < n; j++) {
        prox[j] = std::max(0.0, x[j] - theta);
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
