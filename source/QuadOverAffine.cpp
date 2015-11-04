/* 
 * File:   QuadOverAffine.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on July 24, 2015, 4:55 PM
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

#include "QuadOverAffine.h"
#include "LDLFactorization.h"

void checkConstructorArguments(const Matrix& Q, const Matrix& q, const Matrix& A, const Matrix& b);

QuadOverAffine::QuadOverAffine() : Function() {
    m_A = NULL;
    m_F = NULL;
    m_Fsolver = NULL;
    m_Q = NULL;
    m_q = NULL;
    m_b = NULL;
    m_sigma = NULL;
}

QuadOverAffine::~QuadOverAffine() {
    if (m_Fsolver != NULL) {
        delete m_Fsolver;
    }
    if (m_F != NULL) {
        delete m_F;
    }
    if (m_sigma != NULL) {
        delete m_sigma;
    }
}

void checkConstructorArguments(const Matrix& Q, const Matrix& q, const Matrix& A, const Matrix& b) {
    if (Q.getNrows() != Q.getNcols()) {
        throw std::invalid_argument("Matrix Q is not square");
    }
    if (!q.isColumnVector()) {
        throw std::invalid_argument("q is not a column vector");
    }
    if (!b.isColumnVector()) {
        throw std::invalid_argument("b is not a column vector");
    }
    if (Q.getNcols() != q.getNrows()) {
        throw std::invalid_argument("Q and q have incompatible dimensions");
    }
    if (Q.getNcols() != A.getNcols()) {
        throw std::invalid_argument("Q and A have incompatible dimensions");
    }
    if (A.getNrows() != b.getNrows()) { // A and b must have the same number of rows
        throw std::invalid_argument("A and b have incompatible dimensions");
    }
}

QuadOverAffine::QuadOverAffine(Matrix& Q, Matrix& q, Matrix& A, Matrix& b) {
    checkConstructorArguments(Q, q, A, b);
        
    m_F = NULL;
    m_Fsolver = NULL;
    m_sigma = NULL;
    
    this->m_Q = &Q;
    this->m_q = &q;
    this->m_A = &A;
    this->m_b = &b;    
    size_t n = Q.getNrows();
    size_t s = A.getNrows();
    size_t nF = n + s;
    if (Q.getType() == Matrix::MATRIX_DENSE) {
        m_F = new Matrix(nF, nF, Matrix::MATRIX_DENSE);
        /*
         * F = [Q  * ; *  *]
         */
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                m_F->set(i, j, Q.get(i, j));
            }
        }
        /*
         * F = [Q  A' ; A  *]
         */
        for (size_t i = 0; i < s; i++) {
            for (size_t j = 0; j < n; j++) {
                m_F->set(i + n, j, A.get(i, j));
                m_F->set(j, i + n, A.get(i, j));
            }
        }
    }
    if (m_F != NULL) {
        m_Fsolver = new LDLFactorization(*m_F);
        int status = m_Fsolver -> factorize();
        if (ForBESUtils::STATUS_OK != status) {
            throw std::invalid_argument("LDL factorization failed for matrix F = [Q A'; A 0] (dense) - invalid arguments Q and A");
        }
    }
    m_sigma = new Matrix(nF, 1, Matrix::MATRIX_DENSE);
    for (size_t i = 0; i < s; i++) {
        m_sigma->set(i + n, 0, b.get(i, 0));
    }
}

int QuadOverAffine::callConj(const Matrix& y, double& f_star, Matrix& grad) {
    /* Update sigma(x) */
    for (size_t i = 0; i < m_Q->getNrows(); i++) {
        m_sigma->set(i, 0, y.get(i, 0) - m_q->get(i, 0));
    }
    /* Solve F*grad = sigma */
    int status = m_Fsolver->solve(*m_sigma, grad);
    /* Take the first n elements of grad */
    grad.reshape(m_Q->getNrows(), 1);
    /* f_star = grad' * Q * grad / 2.0 */
    f_star = m_Q->quad(grad);
    for (size_t i = 0; i < grad.getNrows(); i++) { /* Dot product */
        f_star += grad.get(i, 0) * (m_q->get(i, 0) - y.get(i, 0));
    }
    f_star = -f_star;
    return status;
}


