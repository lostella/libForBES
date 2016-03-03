/* 
 * File:   QuadraticLossOverAffine.h
 * Author: chung
 *
 * Created on November 3, 2015, 3:58 PM
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

#ifndef QUADRATICLOSSOVERAFFINE_H
#define	QUADRATICLOSSOVERAFFINE_H

#include "Function.h"
#include "FactoredSolver.h"
#include "LDLFactorization.h"
#include "S_LDLFactorization.h"
#include <math.h>

#define __QUADLOSS_AFFINE_EPSILON 1e-6

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
 * and we perform once and store the LDL-factorization of \f$FF^{\top} + \epsilon I\f$
 * for a small \f$\epsilon>0\f$, that is
 * 
 * \f[
 *  FF^{\top} + \epsilon I = LDL^{\top}.
 * \f]
 * 
 * For a given \f$y\f$ we define the vector \f$\sigma=\sigma(y)\f$ as
 * 
 * \f[
 *  \sigma_i = \frac{y_i}{w_i} + p_i,
 * \f]
 * 
 * and let \f$q=q(y)\f$ be the solution of the linear system
 * 
 * \f[
 *  (FF^{\top} + \epsilon I)q = A\sigma - b.
 * \f]
 * 
 * We solve this system using the above LDL-factorization.
 * 
 * Now the gradient of the conjugate function is given by
 * 
 * \f[
 *  \nabla f^*(y)_i = \sigma_i  - \frac{(A^{\top}q)_i}{w_i}.
 * \f]
 * 
 * and \f$f^*(y)\f$ is given by
 * 
 * \f[
 * f^*(y) = y^{\top}\nabla f^*(y) - \frac{1}{2}\bar{g}^{\top}\mathrm{diag}(w_i)\bar{g},
 * \f]
 * 
 * with \f$\bar{g} = \nabla f^*(y) - p\f$.
 * 
 * \sa S_LDLFactorization
 * \sa FactoredSolver
 * \sa QuadraticLoss
 * 
 * 
 */
class QuadraticLossOverAffine : public Function {
public:
    
    using Function::callConj;
    
    /**
     * Creates a new instance of QuadraticLossOverAffine providing the pair
     * \f$(A,b)\f$ which defines the affine space \f$\mathcal{Z} = 
     * \{z\in\mathbb{R}^n \mid A z = b\}\f$ and the function parameters \f$w\f$
     * and \f$p\f$.
     * 
     * @param A %Matrix A in the definition of the affine space \f$\mathcal{Z}\f$
     * @param b Vector b in the definition of the affine space \f$\mathcal{Z}\f$
     * @param w %Function parameter w
     * @param p %Function parameter p
     * 
     * \exception std::invalid_argument if the parameters have incompatible
     * dimensions.
     * 
     * \exception std::invalid_argument in case the matrix \f$F = AA^{\top} + epsilon I\f$
     * cannot be factorized.
     */
    QuadraticLossOverAffine(Matrix& A, Matrix& b, Matrix& w, Matrix& p);
        
    /**
     * Destructor.
     */
    virtual ~QuadraticLossOverAffine();
    

    virtual int callConj(Matrix& y, double& f_star, Matrix& grad);
    
    virtual int callConj(Matrix& y, double& f_star);

    virtual FunctionOntologicalClass category();



private:

    Matrix * m_A;
    Matrix * m_b;
    Matrix * m_w;
    Matrix * m_p;
    Matrix * m_F;
    FactoredSolver * m_solver;
    
};

#endif	/* QUADRATICLOSSOVERAFFINE_H */

