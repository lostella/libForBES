/* 
 * File:   Norm.h
 * Author: chung
 *
 * Created on October 30, 2015, 5:32 PM
 */

#ifndef NORM_H
#define	NORM_H

#include "Function.h"
#include <math.h>

class Norm : public Function {
public:

    /**
     * Define the dual norm which is defined as
     * \f[
     *  \|x\|_* = \sup_{\|z\|\leq 1} \langle z, x\rangle.
     * \f]
     * 
     * Knowing the dual norm allows us to define the conjugate function of a norm, that is
     * for \f$f(x)=\|x\|\f$, its conjugate is given by
     * 
     * \f[
     *  f^*(y) = 
     *  \begin{cases}
     *      0,&\text{if } \|y\|_* \leq 1,\\
     *      \infty,&\text{otherwise}
     *  \end{cases}
     * \f]
     * 
     * The conjugate function of \f$f(x)=\|x\|\f$ is the indicator function of the
     * dual ball \f$\mathcal{B}_* = \{y\mid \|y\|_* \leq 1 \}\f$.
     * 
     * @param x point in \f$\mathbb{R}^n\f$ where the dual norm should be computed
     * @param norm the value of the dual norm at x, that is \f$\|x\|_*\f$
     * @return 
     * status code which is equal to <code>STATUS_OK</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     */
    virtual int dualNorm(const Matrix &x, double &norm);

    virtual int callConj(const Matrix& x, double& f_star);

protected:
    Norm();

    virtual ~Norm();
private:


};

#endif	/* NORM_H */

