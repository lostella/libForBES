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

    int m_status;
    int m_flag_evalFBE;
    int m_flag_gradFBE;

    FBProblem & m_prob;
    Matrix * m_x;

       
    /* Internal storage for computing proximal-gradient steps */
    Matrix * m_y; /**< x-gamma nabla f(x) */
    Matrix * m_z; /**< prox_{gamma g} (y) */
    Matrix * m_FPRx; /**< x-z */
    Matrix * m_res1x;  /**< L1*x + d1 */
    Matrix * m_gradf1x; /**< nabla f_1 (res1x)*/
    Matrix * m_res2x;
    Matrix * m_gradf2x;
    Matrix * m_gradfx; /**< f(x) = f1(L1*x+d1) + f2(L2*x+d2); gradfx = nabla f(x) */
    Matrix * m_gradFBEx; /** gradient of the FB envelope */    
    double m_f1x; /**< f1(res1x) */
    double m_f2x;
    double m_linx; /**< l'*x */
    double m_fx; /**< f(x) */
    double m_gz; 
    double m_gamma;
    double m_FBEx; /**< FBE(x) */
    double m_sqnormFPRx; /**< ||x-z||^2 */

    /**
     * Evaluates \f$f(x)\f$ and updates the internal status of FBCache.
     *
     * @return status code (see ForBESUtils).
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
     * (using set_point) then the status is reset to \c STATUS_NONE; if \c gamma
     * instead is changed, the status is reset to \c STATUS_EVALF. In fact, the
     * value of f is independent of gamma, and is not to be recomputed.
     * 
     * @param status a status code (see static private const members)
     */
    void reset(int status);

public:
    
    /**
     * 
     * 
     * The internal state of FBCache is stale, or it has not been updated.
     */
    static const int STATUS_NONE = 0; 
    /**
     * 
     * 
     * The following are available:
     * * The <em>residuals</em> \f$ r_1(x)=L_1 x + d_1 \f$ and \f$ r_2(x)=L_2 x + d_2 \f$
     * * The product \f$ \langle l, x \rangle \f$
     * * The values \f$ f_1(r_1(x)) \f$ and \f$ f_2(r_2(x)) \f$
     * * The value of \f$ f(x) = f_1(r_1(x)) + f_2(r_2(x))\f$
     * * The gradient of \f$f_1\f$ at the corresponding residual, \f$ \nabla f_1( r_1(x) ) \f$
     * 
     */
    static const int STATUS_EVALF = 1;
    
    /**
     * 
     * 
     * The data which correspond to #STATUS_EVALF are available and, additionally,
     * the following are available
     * * The gradient of \f$f_2\f$ at the corresponding residual, \f$ \nabla f_2( r_2(x) ) \f$
     * * The gradient of \f$f\f$ at \f$x\f$, that is \f$\nabla f(x)\f$
     * * The value of \f$y = x - \gamma \nabla f(x)\f$
     */
    static const int STATUS_FORWARD = 2;
    /**
     * 
     * Everything has been computed and is available, that is, all data which 
     * correspond to #STATUS_FORWARD (and #STATUS_EVALF) and, additionally,
     * * The proximal \f$z=\mathrm{prox}_{\gamma g}(y)\f$ where \f$y = x - \gamma \nabla f(x)\f$     
     * * The value \f$g(x)\f$
     * * The square norm \f$\|x-z\|^2\f$
     */
    static const int STATUS_FORWARDBACKWARD = 3; /** everything */
    
    
    /**
     * Initialize an FBCache object
     *
     * @param p reference to FBProblem
     * @param x reference to Matrix, the point at which to perform the operations
     * @param gamma the initial stepsize parameter for the operations
     */
    FBCache(FBProblem & p, Matrix & x, double gamma);

    /**
     * Default destructor.
     */
    virtual ~FBCache();

    /**
     * Sets the point at which the FBCache object refers
     *
     * @param x new point at which to evaluate the steps
     */
    void set_point(Matrix& x);

    /**
     * Gets (a pointer to) the point \f$x\f$ to which the FBCache object refer
     *
     * @return a pointer to Matrix containing the handled point
     */
    Matrix * get_point();
    
    /**
     * Gets the result of the forward (gradient) step, with step-size \f$\gamma\f$, at \f$x\f$
     *
     * @param gamma stepsize parameter
     * @return a pointer to Matrix containing the forward step
     */
    Matrix * get_forward_step(double gamma);

    /**
     * Gets the result of the forward-backward (proximal-gradient) with step-size 
     * \f$\gamma\f$ at \f$x\f$, that is \f$z=\mathrm{prox}_{\gamma g}(y)\f$ where \f$y = x - \gamma \nabla f(x)\f$ 
     *
     * @param gamma step-size parameter
     * @return a pointer to Matrix containing the forward-backward step
     */
    Matrix * get_forward_backward_step(double gamma);

    /**
     * Gets the fixed-point residual at \f$x\f$ with parameter \f$\gamma\f$
     *
     * @return a pointer to Matrix containing the fixed-point residual
     */
    Matrix * get_fpr();

    /**
     * Gets the norm of the fixed point residual (FPR) which is defined as
     * \f[
     *   R_{\gamma}(x) = x - \mathrm{prox}_{\gamma g}(x - \gamma \nabla f(x)),
     * \f]
     * where \f$f\f$ is assume to to differentiable and \f$g\f$ is such that it
     * is possible to compute its proximal.
     *
     * @return Euclidean norm of the fixed-point residual, that is \f$R_{\gamma}(x)\f$
     */
    double get_norm_fpr();

    /**
     * Gets the value of \f$f\f$ at \f$x\f$
     *
     * @return \f$f(x)\f$
     */
    double get_eval_f();

    /**
     * Gets the value of the forward-backward envelope function (FBE) 
     * at the current value of \f$x\f$ with a given parameter \f$\gamma\f$
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
     * Erases the internal status of the cache, i.e., sets its status to
     * STATUS_NONE. This means that any getter will require to recompute all 
     * steps.
     */
    void reset();
};

#endif /* FBCACHE_H */
