/* 
 * File:   Quadratic.cpp
 * Author: Chung
 * 
 * Created on July 9, 2015, 3:36 AM
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

#include "Quadratic.h"

using namespace std;

Quadratic::Quadratic() {
    is_Q_eye = true;
    is_q_zero = true;
    L = NULL;
    Q = NULL;
    q = NULL;
}

Quadratic::Quadratic(Matrix& QQ) {
    Q = &QQ;
    is_Q_eye = false;
    is_q_zero = true;
    L = NULL;
    q = NULL;
}

Quadratic::Quadratic(Matrix& QQ, Matrix& qq) {
    q = &qq;
    Q = &QQ;
    L = NULL;
    is_Q_eye = false;
    is_q_zero = false;
}

Quadratic::Quadratic(const Quadratic& orig) {
    // copy-constructor
}

Quadratic::~Quadratic() {
    if (L != NULL) {
        delete L;
    }
}

void Quadratic::setQ(Matrix& Q) {
    is_Q_eye = false;
    this->Q = &Q;
    this->L = NULL;
}

void Quadratic::setq(Matrix& q) {
    is_q_zero = false;
    this->q = &q;
}

int Quadratic::call(Matrix& x, double& f) {
    if (!is_Q_eye) {
        if (is_q_zero) {
            f = Q->quad(x);
        } else {
            f = Q->quad(x, *q);
        }
    } else {
        f = (x * x).get(0, 0);
        if (!is_q_zero) {
            f += ((*q) * x).get(0, 0);
        }
    }
    return ForBESUtils::STATUS_OK;
}

int Quadratic::category() {
    return CAT_QUADRATIC;
}

int Quadratic::callConj(const Matrix& y, double& f_star) {
    Matrix g;
    int status = callConj(y, f_star, g);
    return status;
}

int Quadratic::callConj(const Matrix& y, double& f_star, Matrix& g) {
    Matrix z = (is_q_zero || q == NULL) ? y : y - *q; // z = y    

    if (is_Q_eye || Q == NULL) {
        g = z;
        f_star = (z * z).get(0, 0);
        return ForBESUtils::STATUS_OK;
    }

    if (Q != NULL && Matrix::MATRIX_DIAGONAL == Q->getType()) {
        /* Q is diagonal */
        g = z;
        f_star = 0.0;
        for (size_t i = 0; i < z.getNrows(); i++) {
            g.set(i, 0, g.get(i, 0) / Q->get(i, i));
            f_star += z.get(i, 0) * g.get(i, 0);
        }
        return ForBESUtils::STATUS_OK;
    }

    if (L == NULL) {
        L = new Matrix();
        int status = Q->cholesky(*L);
        if (0 != status) {
            return ForBESUtils::STATUS_NUMERICAL_PROBLEMS;
        }
    }

    L->solveCholeskySystem(g, z); // g = Q \ z 
    f_star = (z * g).get(0, 0); // fstar = z' *g 
    return ForBESUtils::STATUS_OK;
}

int Quadratic::computeGradient(Matrix& x, Matrix& grad) {
    if (is_Q_eye) {
        grad = x;
    } else {
        grad = (*Q) * x;
    }
    if (!is_q_zero) {
        grad += *q;
    }
    return ForBESUtils::STATUS_OK;
}
