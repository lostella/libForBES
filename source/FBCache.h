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
 * \brief Low-level forward-backward operations management class
 * \version version 0.0
 * \ingroup FBSolver-group
 * 
 * FBCache handles forward-backward operations related to an FBProblem
 * In particular, objects of the FBCache class are initialized given an
 * FBProblem p, a point x (of class Matrix) and a scalar parameter gamma.
 * It allows to evaluate the proximal-gradient operation starting from x,
 * and the value of the forward-backward envelope function (FBE) associated
 * with p.
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

    /**
     * Evaluates f(x) and updates the internal status
     *
     * @return Status code, see ForBESUtils.
     */
    int update_eval_f();

    /**
     * Evaluates the forward (gradient) step at x with parameter gamma,
     * and updates the internal status.
     *
     * @return Status code, see ForBESUtils.
     */
    int update_forward_step(double gamma);

    /**
     * Evaluates the forward-backward (or proximal-gradient) step at x with parameter gamma,
     * and updates the internal status.
     *
     * @return Status code, see ForBESUtils.
     */
    int update_forward_backward_step(double gamma);

    /**
     * Evaluates the FBE at x with parameter gamma,
     * and updates the internal status.
     *
     * @return Status code, see ForBESUtils.
     */
    int update_eval_FBE(double gamma);

    /**
     * Evaluates the gradient of the FBE at x with parameter gamma,
     * and updates the internal status.
     *
     * @return Status code, see ForBESUtils.
     */
    int update_grad_FBE(double gamma);

    /**
     * Used to set the internal status of the object at a specific value.
     * For example, if the point at which to evaluate operations is changed
     * (using set_point) then the status is reset to STATUS_NONE; if gamma
     * instead is changed, the status is reset to STATUS_EVALF. In fact, the
     * value of f is independent of gamma, and is not to be recomputed.
     * 
     * @param status a status code (see static private const members)
     */
    void reset(int status);

public:
    /**
     * Initialize an FBCache object
     *
     * @param p reference to FBProblem
     * @param x reference to Matrix, the point at which to perform the operations
     * @param gamma the initial stepsize parameter for the operations
     */
    FBCache(FBProblem & p, Matrix & x, double gamma);

    virtual ~FBCache();
    
    /**
     * Sets the point at which the FBCache object refers
     *
     * @param x new point at which to evaluate the steps
     */
    void set_point(Matrix& x);

    /**
     * Gets (a pointer to) the point at which the FBCache object refers
     *
     * @return a pointer to Matrix containing the handled point
     */
    Matrix * get_point();
    
    /**
     * Gets the result of the forward (gradient) step, with stepsize gamma, at x
     *
     * @param gamma stepsize parameter
     * @return a pointer to Matrix containing the forward step
     */
    Matrix * get_forward_step(double gamma);

    /**
     * Gets the result of the forward-backward (proximal-gradient) with stepsize gamma step at x
     *
     * @param gamma stepsize parameter
     * @return a pointer to Matrix containing the forward-backward step
     */
    Matrix * get_forward_backward_step(double gamma);
    
    /**
     * Gets the fixed-point residual at x with parameter gamma
     *
     * @return a pointer to Matrix containing the fixed-point residual
     */
    Matrix * get_fpr();

    /**
     * Gets the norm of the fixed point residual
     *
     * @return Euclidean norm of the fixed-point residual
     */
    double get_norm_fpr();

    /**
     * Gets the value of f at x
     *
     * @return f(x)
     */
    double get_eval_f();

    /**
     * Gets the value of the FBE at x with a given parameter gamma
     *
     * @param gamma stepsize parameter
     * @return value of the FBE at x with parameter gamma
     */
    double get_eval_FBE(double gamma);

    /**
     * Gets the gradient of the FBE at x with a given parameter gamma
     *
     * @param gamma stepsize parameter
     * @return pointer to Matrix containing the gradient of the FBE
     */
    Matrix * get_grad_FBE(double gamma);

    /**
     * Erases the internal status of the cache, i.e., set its status to
     * STATUS_NONE. This means that any get_<something> call will require
     * recomputing all the steps.
     */
    void reset();
};

#endif /* FBCACHE_H */
