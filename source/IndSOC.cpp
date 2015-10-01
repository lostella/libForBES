/*
 * File:   IndSOC.cpp
 * Author: Lorenzo Stella
 *
 * Created on September 21, 2015, 11:08 AM
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
 
#include <algorithm>
#include <math.h>

#include "IndSOC.h"

IndSOC::IndSOC(int n) : Function() {
    this->n = n;
}

IndSOC::~IndSOC() {
}

int IndSOC::category() {
    return Function::CAT_INDICATOR;
}

int IndSOC::call(Matrix& x, double& f) {
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a vector");
    }
    bool isInside = true;
    size_t i = 0;
    size_t n = x.getNrows();
    double squaredNorm = 0, xi;
    double t = x.get(n-1, 0);

    while (i < n-1 && isInside) {
        xi = x.get(i, 0);
        squaredNorm += xi*xi;
        i++;
    }

    f = (t >= sqrt(squaredNorm)) ? 0.0 : INFINITY;
    return ForBESUtils::STATUS_OK;
}

int IndSOC::callProx(const Matrix& x, double gamma, Matrix& prox) {
    double f_at_prox;
    return callProx(x, gamma, prox, f_at_prox);
}

int IndSOC::callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    f_at_prox = 0.0;
    assert(x.isColumnVector());

    size_t i = 0;
    size_t n = x.getNrows();
    double norm = 0, xi, scal;
    double t = x.get(n-1, 0);
//    double abst = t > 0 ? t : -t;

    while (i < x.getNrows()-1) {
        xi = x.get(i, 0);
        norm += xi*xi;
        i++;
    }
    
    norm = sqrt(norm);
    
    if (t > norm) {
        /* return the original point */
        /* is there a faster way to do this? for sure */
        for (i = 0; i < n; i++) {
            prox.set(i, 0, x.get(i, 0));
        }
    } else if (t < -norm) {
        /* return the zero vector */
        /* is there a faster way to do this? for sure */
        for (i = 0; i < n; i++) {
            prox.set(i, 0, 0);
        }
    } else {
        /* perform actual projection here */
        scal = (1+t/norm)/2;
        for (i = 0; i < n-1; i++) {
            xi = x.get(i, 0);
            prox.set(i, 0, scal*xi);
        }
        prox.set(n-1, 0, (norm + t)/2);
    }
    
    return ForBESUtils::STATUS_OK;
}
