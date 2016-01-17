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
#include "MatrixFactory.h"
#include <vector>
#include <algorithm>
#include <cmath>

IndProbSimplex::IndProbSimplex() {

}

IndProbSimplex::~IndProbSimplex() {
}

int IndProbSimplex::callProx(Matrix& x, double gamma, Matrix& prox) {

    size_t n = x.getNrows();

    std::vector<double> x_hat_vec(x.getData(), x.getData() + n);

    /* x_hat := rev_sort(x_hat) */
    std::sort(x_hat_vec.rbegin(), x_hat_vec.rend());

    double t = 0.0;
    bool flag = true;
    size_t i = 0;
    double val;

    /* skipping negative data */
    const std::vector<double>::iterator lb = std::lower_bound(x_hat_vec.begin(), x_hat_vec.end(), 0.0, std::greater<double>());

    std::vector<double>::iterator it;
    it += (lb - x_hat_vec.begin());

    for (it = x_hat_vec.begin(); flag && it != x_hat_vec.end(); ++it, ++i) {
        val = std::max(0.0, *it);
        t += val;
        flag = val * (i + 1.0) > t - 1.0;
    }
    t -= 1.0 + val;

    const double theta = std::max(0.0, t / (i - 1));

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
