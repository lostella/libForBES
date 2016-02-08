/*
 * File:   FBCache.cpp
 * Author: Lorenzo Stella, Pantelis Sopasakis
 *
 * Created on October 2, 2015
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

#include "FBCache.h"
#include "LinearOperator.h"

// #include <iostream>
#include <cmath>
#include <limits>

void FBCache::reset(int status) {
    if (status < m_status) m_status = status;
    m_flag_evalFBE = 0;
    m_flag_gradFBE = 0;
}

void FBCache::reset() {
    FBCache::reset(FBCache::STATUS_NONE);
}

FBCache::FBCache(FBProblem & p, Matrix & x, double gamma) : m_prob(p), m_x(&x), m_gamma(gamma) {
    reset(FBCache::STATUS_NONE);


    // get dimensions of things
    size_t m_x_rows = m_x->getNrows();
    size_t m_x_cols = m_x->getNcols();
    size_t m_res1_rows, m_res1_cols;
    size_t m_res2_rows, m_res2_cols;

    if (m_prob.d1() != NULL) {
        m_res1_rows = m_prob.d1()->getNrows();
        m_res1_cols = m_prob.d1()->getNcols();
    } else if (m_prob.L1() != NULL) {
        m_res1_rows = m_prob.L1()->dimensionOut().first;
        m_res1_cols = m_prob.L1()->dimensionOut().second;
    } else {
        m_res1_rows = m_x_rows;
        m_res1_cols = m_x_cols;
    }
    if (m_prob.d2() != NULL) {
        m_res2_rows = m_prob.d2()->getNrows();
        m_res2_cols = m_prob.d2()->getNcols();
    } else if (m_prob.L2() != NULL) {
        m_res2_rows = m_prob.L2()->dimensionOut().first;
        m_res2_cols = m_prob.L2()->dimensionOut().second;
    } else {
        m_res2_rows = m_x_rows;
        m_res2_cols = m_x_cols;
    }

    // allocate memory for residuals and gradients (where needed)
    if (m_prob.f1() != NULL) {
        m_res1x = new Matrix(m_res1_rows, m_res1_cols);
        m_gradf1x = new Matrix(m_res1_rows, m_res1_cols);
    } else {
        m_res1x = NULL;
        m_gradf1x = NULL;
    }
    if (m_prob.f2() != NULL) {
        m_res2x = new Matrix(m_res2_rows, m_res2_cols);
        m_gradf2x = new Matrix(m_res2_rows, m_res2_cols);
    } else {
        m_res2x = NULL;
        m_gradf2x = NULL;
    }

    m_gradfx = new Matrix(m_x_rows, m_x_cols);
    m_z = new Matrix(m_x_rows, m_x_cols);
    m_y = new Matrix(m_x_rows, m_x_cols);
    m_FPRx = new Matrix(m_x_rows, m_x_cols);
    m_gradFBEx = new Matrix(m_x_rows, m_x_cols);

    m_FBEx = std::numeric_limits<double>::infinity();
    m_sqnormFPRx = std::numeric_limits<double>::infinity();

    m_f1x = 0.0;
    m_f2x = 0.0;
    m_linx = 0.0;
    m_fx = 0.0;
    m_gz = 0.0;
}

int FBCache::update_eval_f() {
    if (m_status >= FBCache::STATUS_EVALF) {
        return ForBESUtils::STATUS_OK;
    }

    if (m_prob.f1() != NULL) {
        if (m_prob.L1() != NULL) {
            *m_res1x = m_prob.L1()->call(*m_x);
        } else {
            *m_res1x = *m_x;
        }
        if (m_prob.d1() != NULL) {
            *m_res1x += *(m_prob.d1());
        }
        int status = m_prob.f1()->call(*m_res1x, m_f1x, *m_gradf1x);
        if (ForBESUtils::STATUS_OK != status) {
            return status;
        }
    }

    if (m_prob.f2() != NULL) {
        if (m_prob.L2() != NULL) {
            *m_res2x = m_prob.L2()->call(*m_x);
        } else {
            *m_res2x = *m_x;
        }
        if (m_prob.d2() != NULL) {
            *m_res2x += *(m_prob.d2());
        }
        int status = m_prob.f2()->call(*m_res2x, m_f2x);
        if (ForBESUtils::STATUS_OK != status) {
            return status;
        }
    }

    if (m_prob.lin() != NULL) {
        m_linx = ((*m_prob.lin()) * (*m_x))[0];
    }

    m_fx = m_f1x + m_f2x + m_linx;
    m_status = FBCache::STATUS_EVALF;

    return ForBESUtils::STATUS_OK;
}

int FBCache::update_forward_step(double gamma) {
    if (gamma != m_gamma) {
        reset(FBCache::STATUS_EVALF);
    }

    if (m_status >= FBCache::STATUS_FORWARD) {
        *m_y = *m_x;
        Matrix::add(*m_y, -gamma, *m_gradfx, 1.0);
        m_gamma = gamma;
        return ForBESUtils::STATUS_OK;
    }
    int status;
    if (m_status < FBCache::STATUS_EVALF) {
        status = update_eval_f();
        if (!ForBESUtils::is_status_ok(status)) {
            return status;
        }
    }

    if (m_prob.f1() != NULL) {
        if (m_prob.L1()) {
            Matrix d_gradfx = m_prob.L1()->callAdjoint(*m_gradf1x);
            *m_gradfx = d_gradfx;
        } else {
            *m_gradfx = *m_gradf1x;
        }
    }

    if (m_prob.f2() != NULL) {
        status = m_prob.f2()->call(*m_x, m_f2x, *m_gradf2x);
        if (!ForBESUtils::is_status_ok(status)) {
            return status;
        }
        if (m_prob.L2() != NULL) {
            Matrix d_gradfx = m_prob.L2()->callAdjoint(*m_gradf2x);
            if (m_prob.f1() != NULL) *m_gradfx += d_gradfx;
            else *m_gradfx = d_gradfx;
        } else {
            if (m_prob.f1() != NULL) *m_gradfx += *m_gradf2x;
            else *m_gradfx = *m_gradf2x;
        }
    }

    if (m_prob.lin()) {
        if (m_prob.f1() != NULL || m_prob.f2() != NULL) {
            *m_gradfx += (*m_prob.lin());
        } else {
            *m_gradfx = *m_prob.lin();
        }
    }

    *m_y = *m_x;
    Matrix::add(*m_y, -gamma, *m_gradfx, 1.0);

    m_gamma = gamma;
    m_status = FBCache::STATUS_FORWARD;

    return ForBESUtils::STATUS_OK;
}

int FBCache::update_forward_backward_step(double gamma) {
    int status;
    if (gamma != m_gamma) {
        reset(FBCache::STATUS_EVALF);
    }
    if (m_status >= FBCache::STATUS_FORWARDBACKWARD) {
        return ForBESUtils::STATUS_OK;
    }
    if (m_status < FBCache::STATUS_FORWARD) {
        status = update_forward_step(gamma);
        if (!ForBESUtils::is_status_ok(status)) {
            return status;
        }
    }
    /*
     * Shouldn't we check whether g() is not NULL
     */
    status = m_prob.g()->callProx(*m_y, gamma, *m_z, m_gz);
    if (!ForBESUtils::is_status_ok(status)) {
        return status;
    }
    *m_FPRx = (*m_x - *m_z);
    m_sqnormFPRx = 0.0;
    for (int i = 0; i < m_FPRx->length(); i++) {
        double c = (*m_FPRx)[i];
        m_sqnormFPRx += c*c;
    }

    m_gamma = gamma;
    m_status = FBCache::STATUS_FORWARDBACKWARD;

    return ForBESUtils::STATUS_OK;
}

int FBCache::update_eval_FBE(double gamma) {
    if (gamma != m_gamma) { /* Use tolerance */
        reset(FBCache::STATUS_EVALF);
    }

    if (m_flag_evalFBE == 1) {
        return ForBESUtils::STATUS_OK;
    }

    if (m_status < FBCache::STATUS_FORWARDBACKWARD) {
        int status = update_forward_backward_step(gamma);
        if (!ForBESUtils::is_status_ok(status)) {
            return status;
        }
    }

    Matrix innprox_mat(1, 1);
    m_FPRx -> transpose();
    Matrix::mult(innprox_mat, 1.0, *m_FPRx, *m_gradfx, 0.0);
    m_FPRx -> transpose();
    double innprod = innprox_mat[0];

    m_FBEx = m_fx + m_gz - innprod + 0.5 / m_gamma*m_sqnormFPRx;
    m_gamma = gamma;
    m_flag_evalFBE = 1;

    return ForBESUtils::STATUS_OK;
}

int FBCache::update_grad_FBE(double gamma) {
    if (gamma != m_gamma) {
        reset(FBCache::STATUS_EVALF);
    }

    if (m_flag_gradFBE == 1) {
        return ForBESUtils::STATUS_OK;
    }

    if (m_status < FBCache::STATUS_FORWARDBACKWARD) {
        int status = update_forward_backward_step(gamma);
        if (!ForBESUtils::is_status_ok(status)) {
            return status;
        }
    }

    *m_gradFBEx = *m_FPRx;

    // gradFBE(x) = (I-gamma*H(x))*FPR(x)/gamma
    // if the smooth term is f(Lx+d) then H(x) = L'*hessf(x)*L
    // so in general the way to compute gradFBE(x) should be:
    //
    //  v1 <- L*FPR(x)
    //  v2 <- H(x)*v1
    //  v3 <- L'*v2
    //  gradFBE(x) <- FPR(x)/gamma - v3
    //
    // and when L is not present (that is L = Identity):
    //
    //  v1 <- H(x)*FPR(x)
    //  gradFBE(x) <- FPR(x)/gamma - v1

    if (m_prob.f1() != NULL) {
        if (m_prob.L1() != NULL) {
            Matrix v1 = m_prob.L1()->call(*m_FPRx);
            Matrix v2 = Matrix(m_prob.L1()->dimensionOut());
            m_prob.f1()->hessianProduct(*m_x, v1, v2);
            Matrix v3 = m_prob.L1()->callAdjoint(v2);
            Matrix::add(*m_gradFBEx, -1.0, v3, 1.0 / gamma);
        } else {
            Matrix v1 = Matrix(m_x->getNrows(), m_x->getNcols());
            m_prob.f1()->hessianProduct(*m_x, *m_FPRx, v1);
            Matrix::add(*m_gradFBEx, -1.0, v1, 1.0 / gamma);
        }
    }

    if (m_prob.f2() != NULL) {
        return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
    }

    m_gamma = gamma;
    m_flag_gradFBE = 1;

    return ForBESUtils::STATUS_OK;
}

void FBCache::set_point(Matrix& x) {
    *m_x = x;
    reset(FBCache::STATUS_NONE);
}

Matrix * FBCache::get_point() {
    return m_x;
}

double FBCache::get_eval_FBE(double gamma) {
    update_eval_FBE(gamma);
    return m_FBEx;
}

Matrix * FBCache::get_grad_FBE(double gamma) {
    update_grad_FBE(gamma);
    return m_gradFBEx;
}

double FBCache::get_eval_f() {
    update_eval_f();
    return m_fx;
}

Matrix* FBCache::get_forward_step(double gamma) {
    update_forward_step(gamma);
    return m_y;
}

Matrix* FBCache::get_forward_backward_step(double gamma) {
    update_forward_backward_step(gamma);
    return m_z;
}

Matrix* FBCache::get_fpr() {
    update_forward_backward_step(m_gamma);
    return m_FPRx;
}

double FBCache::get_norm_fpr() {
    update_forward_backward_step(m_gamma);
    return sqrt(m_sqnormFPRx);
}

FBCache::~FBCache() {
    if (m_z != NULL) {
        delete m_z;
        m_z = NULL;
    }
    if (m_y != NULL) {
        delete m_y;
        m_y = NULL;
    }
    if (m_res2x != NULL) {
        delete m_res2x;
        m_res2x = NULL;
    }
    if (m_gradf2x != NULL) {
        delete m_gradf2x;
        m_gradf2x = NULL;
    }
    if (m_res1x != NULL) {
        delete m_res1x;
        m_res1x = NULL;
    }
    if (m_gradf1x != NULL) {
        delete m_gradf1x;
        m_gradf1x = NULL;
    }
    if (m_gradfx != NULL) {
        delete m_gradfx;
        m_gradfx = NULL;
    }
    if (m_FPRx != NULL) {
        delete m_FPRx;
        m_FPRx = NULL;
    }
    if (m_gradFBEx != NULL) {
        delete m_gradFBEx;
        m_gradFBEx = NULL;
    }
}
