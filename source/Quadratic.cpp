/* 
 * File:   Quadratic.cpp
 * Author: Pantelis Sopasakis
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
    m_is_Q_eye = true;
    m_is_q_zero = true;
    m_solver = NULL;
    m_Q = NULL;
    m_q = NULL;
}

Quadratic::Quadratic(Matrix& QQ) {
    m_Q = &QQ;
    m_is_Q_eye = false;
    m_is_q_zero = true;
    m_solver = NULL;
    m_q = NULL;
}

Quadratic::Quadratic(Matrix& QQ, Matrix& qq) {
    m_q = &qq;
    m_Q = &QQ;
    m_solver = NULL;
    m_is_Q_eye = false;
    m_is_q_zero = false;
}

Quadratic::~Quadratic() {
    if (m_solver != NULL) {
        delete m_solver;
    }
}

void Quadratic::setQ(Matrix& Q) {
    m_is_Q_eye = false;
    this->m_Q = &Q;
    this->m_solver = NULL;
}

void Quadratic::setq(Matrix& q) {
    m_is_q_zero = false;
    this->m_q = &q;
}

int Quadratic::call(Matrix& x, double& f, Matrix& grad) {
    int statusComputeGrad = computeGradient(x, grad); // compute the gradient of f at x (grad)
    if (statusComputeGrad != ForBESUtils::STATUS_OK) {
        return statusComputeGrad;
    }
    // f = (1/2)*(grad+q)'*x
    f = ((m_is_q_zero ? grad : grad + (*m_q)) * x).get(0, 0) / 2;
    return ForBESUtils::STATUS_OK;
}

int Quadratic::call(Matrix& x, double& f) {
    if (!m_is_Q_eye) {
        if (m_is_q_zero) {
            f = m_Q->quad(x);
        } else {
            f = m_Q->quad(x, *m_q);
        }
    } else {
        f = (x * x).get(0, 0);
        if (!m_is_q_zero) {
            f += ((*m_q) * x).get(0, 0);
        }
    }
    return ForBESUtils::STATUS_OK;
}

int Quadratic::callConj(Matrix& y, double& f_star) {
    Matrix g;
    int status = callConj(y, f_star, g);
    return status;
}

int Quadratic::callConj(Matrix& y, double& f_star, Matrix& g) {
    Matrix z = (m_is_q_zero || m_q == NULL) ? y : y - *m_q; // z = y    
    if (m_is_Q_eye || m_Q == NULL) {
        g = z;
        f_star = (z * z).get(0, 0);
        return ForBESUtils::STATUS_OK;
    }
    if (m_Q != NULL && Matrix::MATRIX_DIAGONAL == m_Q->getType()) {
        /* Q is diagonal */
        g = z;
        f_star = 0.0;
        for (size_t i = 0; i < z.getNrows(); i++) {
            g.set(i, 0, g.get(i, 0) / m_Q->get(i, i));
            f_star += z.get(i, 0) * g.get(i, 0);
        }
        return ForBESUtils::STATUS_OK;
    }

    if (m_solver == NULL) {
        m_solver = new CholeskyFactorization(*m_Q);
        int status = m_solver->factorize();
        if (0 != status) {
            return ForBESUtils::STATUS_NUMERICAL_PROBLEMS;
        }
    }

    m_solver->solve(z, g); // Q*g = z   OR  g = Q \ z
    f_star = (z * g).get(0, 0); // fstar = z' *g 
    return ForBESUtils::STATUS_OK;
}

int Quadratic::computeGradient(Matrix& x, Matrix& grad) {
    if (m_is_Q_eye) {
        grad = x;
    } else {
        grad = (*m_Q) * x;
    }
    if (!m_is_q_zero) {
        grad += *m_q;
    }
    return ForBESUtils::STATUS_OK;
}

FunctionOntologicalClass Quadratic::category() {
    return FunctionOntologyRegistry::quadratic();
}

