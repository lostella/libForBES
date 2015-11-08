/* 
 * File:   QuadOverAffine.h
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

#ifndef QUADOVERAFFINE_H
#define	QUADOVERAFFINE_H

#include "Function.h"
#include "FactoredSolver.h"

/**
 * \class QuadOverAffine
 * \brief %Function <code>F(x) = 0.5*x'*Q*x + q'*x + delta(x|Z)</code>
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on July 24, 2015, 4:55 PM
 * 
 * Quadratic-over-affine function.
 * 
 * A quadratic-over-affine function, or, a quadratic function plus the indicator of an affine subspace, 
 * is a function \f$F:\mathbb{R}^n \to \mathbb{R}\cup\{+\infty\}\f$ is the form
 * 
 * \f[
 * F(x) = \frac{1}{2}x'Qx + q'x + \delta(x|E),
 * \f]
 * 
 * where \f$E\f$ is the affine space
 * 
 * \f[
 * E = \{x: Ax = b\}
 * \f]
 * 
 * and \f$A\in\mathbb{R}^{s\times n}\f$, \f$b\in\mathbb{R}^s\f$ and \f$\delta(\cdot|E)\f$ is the indicator function
 * 
 * \f[
 * \delta(x|E) = \begin{cases}
 * 1, \text{ if } x \in E,\\
 * \infty, \text{ otherwise} 
 * \end{cases}
 * \f]
 * 
 * Let us define the matrix
 * 
 * \f[
 * S = \begin{bmatrix}
 * Q & A'\\
 * A & 0
 * \end{bmatrix},
 * \f]
 * 
 * then, the conjugate of <code>F</code> is
 * 
 * \f[
 * S g(x^*) =  \begin{bmatrix}x^*-q\\b\end{bmatrix}
 * \f]
 * and let \f$\gamma(x^*)\f$ be the vector comprising the first n entries of \f$g(x^*)\f$. Then,
 * \f[
 * F^*(x^*) = -\frac{1}{2} \left(\gamma(x^*)'Q\gamma(x^*) + (q-x^*)'\gamma(x^*)\right)
 * \f]
 * 
 * Here is an example of use
 * 
 * \code{.cpp}
 * // First, define matrices Q, q, A and b
 * Matrix Q = ...;
 * Matrix q = ...;
 * Matrix A = ...;
 * Matrix b = ...;
 * Function *F = new QuadOverAffine(Q, q, A, b); // define the function as QuadOverAffine
 * Matrix y = ...;
 * double f_star;
 * Matrix grad;
 * int status = F->callConj(y, f_star, grad);
 * \endcode
 * 
 * \ingroup Functions
 */
class QuadOverAffine : public Function {
public:
    
    using Function::callConj;

    /**
     * Define a new quadratic-over-affine function.
     * 
     * @param Q %Matrix Q
     * @param q Vector q
     * @param A %Matrix A
     * @param b Vector b
     * 
     * \exception std::invalid_argument in case the given parameters have incompatible
     * dimensions or matrix F = [Q A'; A 0] cannot be LDL-decomposed.
     */
    QuadOverAffine(Matrix& Q, Matrix& q, Matrix& A, Matrix& b);

    /**
     * Destructor
     */
    virtual ~QuadOverAffine();

    virtual int callConj(Matrix& y, double& f_star, Matrix& grad);

    virtual FunctionOntologicalClass category();



private:

    QuadOverAffine();

    Matrix *m_Q; /**< Matrix Q (Hessian) */
    Matrix *m_q; /**< Vector q (Linear term) */
    Matrix *m_A; /**< Matrix A */
    Matrix *m_b; /**< Matrix b */

    Matrix *m_F; /**< Matrix <code>F = [Q A'; A 0]</code> */
    Matrix *m_sigma;
    FactoredSolver * m_Fsolver; /**< Factorizer for matrix F */
};

#endif	/* QUADOVERAFFINE_H */

