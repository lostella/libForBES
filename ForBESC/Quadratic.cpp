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

int Quadratic::call( Matrix& x, double& f)  {
    if (!is_Q_eye) {
        if (is_q_zero) {
            f = Q->quad(x);
        } else {
            f = Q->quad(x, *q);
        }
    }
    return STATUS_OK;
}

int Quadratic::category() {
    return CAT_QUADRATIC;
}

int Quadratic::callConj(const Matrix& y, double& f_star) {
    //TODO: Make Cholesky factor (if it doesn't exist)
    if (L == NULL) {
        L = new Matrix();
        int status = Q->cholesky(*L);
        if (0 != status) {
            return STATUS_NUMERICAL_PROBLEMS;
        }
    }
    Matrix z = y - *q;                      // z = y - q
    Matrix g; 
    L->solveCholeskySystem(g, z);           // g = Q \ z 
    f_star = (z * g)[0];                    // fstar = z' *g 
    return STATUS_OK;
}

int Quadratic::callProx(const Matrix& x, double gamma, Matrix& prox, double f_at_prox) {
    return STATUS_UNDEFINED_FUNCTION;
}

int Quadratic::callProx(const Matrix& x, double gamma, Matrix& prox) {
    return STATUS_UNDEFINED_FUNCTION;
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
    return STATUS_OK;
}
