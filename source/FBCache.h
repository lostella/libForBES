/*
 * File:   FBCache.h
 * Author: Lorenzo Stella
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

#ifndef FBCACHE_H
#define	FBCACHE_H

#include "Matrix.h"
#include "FBProblem.h"

/**
 * \class FBCache
 * \brief Low-level cache management system
 * \version version 0.0
 * \ingroup FBSolver-group
 * 
 * 
 */
class FBCache {
private:
	static const int STATUS_NONE = 0;
    static const int STATUS_EVALF = 1;
    static const int STATUS_FORWARD = 2;
    static const int STATUS_FORWARDBACKWARD = 3;
    
    int m_status;
    int m_flag_evalFBE;
    int m_flag_gradFBE;

    FBProblem & m_prob;
    Matrix * m_x;

    /* Reference to problem details */
    Function * m_f1;
    Function * m_f2;
    Function * m_g;
    LinearOperator * m_L1;
    LinearOperator * m_L2;
    Matrix * m_d1;
    Matrix * m_d2;
    Matrix * m_lin;

    /* Vectors dimensions */
    size_t m_x_rows;
    size_t m_x_cols;
    size_t m_res1_rows, m_res1_cols;
    size_t m_res2_rows, m_res2_cols;

    /* Internal storage for computing proximal-gradient steps */
    Matrix * m_y;
    Matrix * m_z;
    Matrix * m_FPRx;
    Matrix * m_res1x;
    Matrix * m_gradf1x;
    Matrix * m_res2x;
    Matrix * m_gradf2x;
    Matrix * m_gradfx;
    Matrix * m_gradFBEx;
    double m_f1x;
    double m_f2x;
    double m_linx;
    double m_fx;
    double m_gz;
    double m_gamma;
    double m_FBEx;
    double m_sqnormFPRx;

    int update_eval_f();
    int update_forward_step(double gamma);
    int update_forward_backward_step(double gamma);

    int update_eval_FBE(double gamma);
    int update_grad_FBE(double gamma);

    void reset(int status);

public:
    FBCache(FBProblem & p, Matrix & x, double gamma);

    virtual ~FBCache();
    
    void set_x(Matrix& x);
    
    Matrix * get_forward_step(double gamma);
    Matrix * get_forward_backward_step(double gamma);

    double get_eval_f();
    double get_eval_FBE(double gamma);
};

#endif /* FBCACHE_H */
