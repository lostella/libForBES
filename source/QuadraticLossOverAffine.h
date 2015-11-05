/* 
 * File:   QuadraticLossOverAffine.h
 * Author: chung
 *
 * Created on November 3, 2015, 3:58 PM
 */

#ifndef QUADRATICLOSSOVERAFFINE_H
#define	QUADRATICLOSSOVERAFFINE_H

#include "Function.h"
#include "FactoredSolver.h"
#include <math.h>


/**
 * \class QuadraticLossOverAffine
 * \brief %Quadratic loss over an affine subspace
 * \version version 0.1
 * \ingroup Functions
 * \date Created on November 3, 2015, 3:58 PM
 * \author Pantelis Sopasakis
 * 
 * We define the quadratic loss over an affine subspace as the following function
 * \f$f:\mathbb{R}^n\to\mathbb{R}\cup \{+\infty\}\f$ given by
 * 
 * \f[
 *   f(x) = \frac{1}{2} \sum_{i=1}^{n}w_i (x_i - p_i)^2 + \delta(x\mid \mathcal{Z}),
 * \f]
 * 
 * where \f$\mathcal{Z}\f$ is an affine subspace of \f$\mathbb{R}^n\f$ given by
 * 
 * \f[
 *  \mathcal{Z} = \{z\in\mathbb{R}^n \mid A z = b\}.
 * \f]
 * 
 * For this function we define its conjugate \f$f^*(y)\f$ and the gradient of
 * its conjugate \f$\nabla f^*(y)\f$.
 * 
 * We first define the matrix \f$F\f$ as
 * 
 * \f[
 * F = A \cdot \mathrm{diag}(1/\sqrt{w_i})_i,
 * \f]
 * 
 * and we perform once and store the LDL-factorization of \f$FF' + \epsilon I\f$
 * for a small \f$\epsilon>0\f$, that is
 * 
 * \f[
 *  FF' + \epsilon I = LDL'.
 * \f]
 * 
 * For a given \f$y\f$ we define the vector \f$\sigma(y)\f$ as
 * 
 * \f[
 *  \sigma_i = \frac{y_i}{w_i} + p_i,
 * \f]
 * 
 * and let \f$q\f$ be the solution of the linear system
 * 
 * \f[
 *  (LDL')q = A\sigma - b.
 * \f]
 * 
 * We solve this system using the above LDL-factorization.
 * 
 * Now the gradient of the conjugate function is given by
 * 
 * \f[
 *  \nabla f^*(y)_i = \sigma_i  - \frac{(A'q)_i}{w_i}.
 * \f]
 * 
 * and \f$f^*(y)\f$ is given by
 * 
 * \f[
 * f^*(y) = y'\nabla f^*(y) - \frac{1}{2}\bar{g}'\mathrm{diag}(w_i)\bar{g},
 * \f]
 * 
 * with \f$\bar{g} = \nabla f^*(y) - p\f$.
 * 
 * \sa LDLFactorization_AAt
 * \sa FactoredSolver
 * \sa QuadraticLoss
 * 
 * 
 */
class QuadraticLossOverAffine : public Function {
public:
    
    QuadraticLossOverAffine(Matrix& A, Matrix& b, Matrix& w, Matrix& p);
    
    
    virtual ~QuadraticLossOverAffine();
private:

    Matrix * m_A;
    Matrix * m_b;
    Matrix * m_w;
    Matrix * m_p;
    FactoredSolver * m_solver;
    
};

#endif	/* QUADRATICLOSSOVERAFFINE_H */

