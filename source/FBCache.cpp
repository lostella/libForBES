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

//#include <iostream>
#include <cmath>

void FBCache::reset(int status) {
    if (status < m_status) m_status = status;
    m_flag_evalFBE = 0;
    m_flag_gradFBE = 0;
}

FBCache::FBCache(FBProblem & p, Matrix & x, double gamma) : m_prob(p), m_x(&x), m_gamma(gamma) {
    reset(FBCache::STATUS_NONE);

    // store pointers to problem and all relevant details
    m_f1 = p.f1();
    m_L1 = p.L1();
    m_d1 = p.d1();
    m_f2 = p.f2();
    m_L2 = p.L2();
    m_d2 = p.d2();
    m_lin = p.lin();
    m_g = p.g();

    // get dimensions of things
    m_x_rows = m_x->getNrows();
    m_x_cols = m_x->getNcols();
    if (m_d1 != NULL) {
        m_res1_rows = m_d1->getNrows();
        m_res1_cols = m_d1->getNcols();
    } else if (m_L1 != NULL) {
        std::pair<size_t, size_t> res1_size = m_L1->dimensionOut();
        m_res1_rows = res1_size.first;
        m_res1_cols = res1_size.second;
    } else {
        m_res1_rows = m_x_rows;
        m_res1_cols = m_x_cols;
    }
    if (m_d2 != NULL) {
        m_res2_rows = m_d2->getNrows();
        m_res2_cols = m_d2->getNcols();
    } else if (m_L2 != NULL) {
        std::pair<size_t, size_t> res2_size = m_L2->dimensionOut();
        m_res2_rows = res2_size.first;
        m_res2_cols = res2_size.second;
    } else {
        m_res2_rows = m_x_rows;
        m_res2_cols = m_x_cols;
    }

    // allocate memory for residuals and gradients (where needed)
    if (m_f1 != NULL) {
        m_res1x = new Matrix(m_res1_rows, m_res1_cols);
        m_gradf1x = new Matrix(m_res1_rows, m_res1_cols);
    } else {
        m_res1x = NULL;
        m_gradf1x = NULL;
    }
    if (m_f2 != NULL) {
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

    if (m_f1 != NULL) {
        //    	cout << "there's f1" << endl;
        if (m_L1 != NULL) {
            //        	cout << "there's L1" << endl;
            *m_res1x = m_L1->call(*m_x);
        } else {
            *m_res1x = *m_x;
        }
        if (m_d1 != NULL) {
            //        	cout << "there's d1" << endl;
            *m_res1x += *m_d1;
        }
        int status = m_f1->call(*m_res1x, m_f1x, *m_gradf1x);
        if (ForBESUtils::STATUS_OK != status) {
            cout << "ERROR: " << status << endl;
            return status;
        }
    }

    if (m_f2 != NULL) {
        if (m_L2 != NULL) *m_res2x = m_L2->call(*m_x);
        else *m_res2x = *m_x;
        if (m_d2 != NULL) *m_res2x += *m_d2;
        int status = m_f2->call(*m_res2x, m_f2x);
        if (ForBESUtils::STATUS_OK != status) {
            return status;
        }
    }

    if (m_lin != NULL) {
        m_linx = ((*m_lin) * (*m_x))[0];
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

    if (m_status < FBCache::STATUS_EVALF) {
        int status = update_eval_f();
    }

    if (m_f1 != NULL) {
        if (m_L1) {
            Matrix d_gradfx = m_L1->callAdjoint(*m_gradf1x);
            *m_gradfx = d_gradfx;
        } else {
            *m_gradfx = *m_gradf1x;
        }
    }

    if (m_f2 != NULL) {
        m_f2->call(*m_x, m_f2x, *m_gradf2x);
        if (m_L2) {
            Matrix d_gradfx = m_L2->callAdjoint(*m_gradf2x);
            if (m_f1 != NULL) *m_gradfx += d_gradfx;
            else *m_gradfx = d_gradfx;
        } else {
            if (m_f1 != NULL) *m_gradfx += *m_gradf2x;
            else *m_gradfx = *m_gradf2x;
        }
    }

    if (m_lin) {
        if (m_f1 != NULL || m_f2 != NULL) *m_gradfx += (*m_lin);
        else *m_gradfx = *m_lin;
    }

    *m_y = *m_x;
    Matrix::add(*m_y, -gamma, *m_gradfx, 1.0);

    m_gamma = gamma;
    m_status = FBCache::STATUS_FORWARD;

    return ForBESUtils::STATUS_OK;
}

int FBCache::update_forward_backward_step(double gamma) {
    if (gamma != m_gamma) {
        reset(FBCache::STATUS_EVALF);
    }

    if (m_status >= FBCache::STATUS_FORWARDBACKWARD) {
        return ForBESUtils::STATUS_OK;
    }

    if (m_status < FBCache::STATUS_FORWARD) {
        int status = update_forward_step(gamma);
    }
    double c;
    int status = m_g->callProx(*m_y, gamma, *m_z, m_gz);
    *m_FPRx = (*m_x - *m_z);
    m_sqnormFPRx = 0;
    for (int i = 0; i < m_FPRx->length(); i++) {
        c = (*m_FPRx)[i];
        m_sqnormFPRx += c*c;
    }

    m_gamma = gamma;
    m_status = FBCache::STATUS_FORWARDBACKWARD;

    return ForBESUtils::STATUS_OK;
}

int FBCache::update_eval_FBE(double gamma) {
    if (gamma != m_gamma) {
        reset(FBCache::STATUS_EVALF);
    }

    if (m_flag_evalFBE == 1) {
        return ForBESUtils::STATUS_OK;
    }

    if (m_status < FBCache::STATUS_FORWARDBACKWARD) {
        int status = update_forward_backward_step(gamma);
    }

    Matrix innprox_mat(1,1);
    m_FPRx -> transpose();
    Matrix::mult(innprox_mat, 1.0, *m_FPRx, *m_gradfx, 0.0);
    m_FPRx -> transpose();
    double innprod = innprox_mat.get(0,0);


    m_FBEx = m_fx + m_gz - innprod + 0.5 / m_gamma*m_sqnormFPRx;

    m_gamma = gamma;
    m_flag_evalFBE = 1;

    return ForBESUtils::STATUS_OK;
}

int FBCache::update_grad_FBE(double gamma) {
    if (gamma != m_gamma) {
        reset(FBCache::STATUS_EVALF);
    }

    if (m_flag_evalFBE == 1) {
        return ForBESUtils::STATUS_OK;
    }

    if (m_status < FBCache::STATUS_FORWARDBACKWARD) {
        int status = update_forward_backward_step(gamma);
    }

    /* TODO: fill in here */

    m_gamma = gamma;
    m_flag_gradFBE = 1;

    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

void FBCache::set_point(Matrix& x) {
    *m_x = x;
    reset(FBCache::STATUS_NONE);
}

double FBCache::get_eval_FBE(double gamma) {
    int status = update_eval_FBE(gamma);
    return m_FBEx;
}

double FBCache::get_eval_f() {
    int status = update_eval_f();
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
}
