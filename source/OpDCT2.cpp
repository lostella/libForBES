/* 
 * File:   OpDCT2.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on September 15, 2015, 3:36 PM
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

#include "OpDCT2.h"
#include <cmath>
#include <iostream>


void update_y_helper_n_even(Matrix& y, double alpha, Matrix& x, double gamma, size_t n);
void update_y_helper_n_odd(Matrix& y, double alpha, Matrix& x, double gamma, size_t n);
double power_of_minus_one(size_t k);

OpDCT2::OpDCT2() : LinearOperator(), m_dimension(_EMPTY_OP_DIM) {

}

OpDCT2::OpDCT2(size_t n) : m_dimension(_VECTOR_OP_DIM(n)) {
}

OpDCT2::~OpDCT2() {
}

static const double FOO[4] = {1.0, 0.0, -1.0, 0.0};

inline double power_of_minus_one(size_t k) {
    return (k % 2 == 0) ? 1.0 : -1.0;
}

void update_y_helper_n_even(Matrix& y, double alpha, Matrix& x, double gamma, size_t n) {
    size_t nu = n / 2;
    bool apply_trick = (n >= 8) && (n % 4 == 0);
    for (size_t k = 1; k < n; k++) {
        double yk = 0.0;
        if (apply_trick && k % 2 == 0) {
            for (size_t i = 0; i < nu / 2; i++) {
                double aik;
                size_t mu = k / 2;
                aik = std::cos(M_PI * (static_cast<double> (i) + 0.5) * static_cast<double> (k) / static_cast<double> (n));
                yk += (
                        x.get(i, 0)
                        + power_of_minus_one(mu) * (x.get(nu - i - 1, 0) + x.get(nu + i, 0))
                        + x.get(n - i - 1, 0)
                        ) * aik;
            }
        } else {
            for (size_t i = 0; i < nu; i++) {
                double aik;
                aik = std::cos(M_PI * (static_cast<double> (i) + 0.5) * static_cast<double> (k) / static_cast<double> (n));
                if (k % 2 == 1) {
                    yk += (x.get(i, 0) - x.get(n - i - 1, 0)) * aik;
                } else {
                    yk += (x.get(i, 0) + x.get(n - i - 1, 0)) * aik;
                }
            }
        }
        y.set(k, 0, gamma * y.get(k, 0) + alpha * yk);
    }
}

void update_y_helper_n_odd(Matrix& y, double alpha, Matrix& x, double gamma, size_t n) {
    for (size_t k = 1; k < n; k++) {
        double yk = 0.0;
        for (size_t i = 0; i < n / 2; i++) {
            double aik;
            aik = std::cos(M_PI * (static_cast<double> (i) + 0.5) * static_cast<double> (k) / static_cast<double> (n));
            yk += (x.get(i, 0) + power_of_minus_one(k) * x.get(n - 1 - i, 0)) * aik;
        }
        yk += FOO[k % 4] * x.get(n / 2, 0);
        y.set(k, 0, gamma * y.get(k, 0) + alpha * yk);
    }
}

int OpDCT2::call(Matrix& y, double alpha, Matrix& x, double gamma) {
    size_t n = x.getNrows();
    double y0 = gamma * y.get(0, 0);
    for (size_t i = 0; i < n; i++) {
        y0 += alpha * x.get(i, 0);
    }
    y.set(0, 0, y0);
    if (n % 2 == 0) { // if n is even
        update_y_helper_n_even(y, alpha, x, gamma, n);
    } else {
        update_y_helper_n_odd(y, alpha, x, gamma, n);
    }
    return ForBESUtils::STATUS_OK;
}

int OpDCT2::callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma) {
    size_t n = x.getNrows();
    for (size_t k = 0; k < n; k++) {
        double v = 0.0;
        for (size_t i = 0; i < n; i++) {
            double aki;
            aki = std::cos(M_PI * (static_cast<double> (k) + 0.5) * static_cast<double> (i) / static_cast<double> (n));
            v += x.get(i, 0) * aki;
        }
        y.set(k, 0, gamma * y.get(k, 0) + alpha * v);
    }
    return ForBESUtils::STATUS_OK;
}

std::pair<size_t, size_t> OpDCT2::dimensionIn() {
    return m_dimension;
}

std::pair<size_t, size_t> OpDCT2::dimensionOut() {
    return m_dimension;
}

bool OpDCT2::isSelfAdjoint() {
    return false;
}
