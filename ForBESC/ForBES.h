/* 
 * File:   ForBES.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 22, 2015, 12:25 PM
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
 * 
 */

#ifndef FORBES_H
#define	FORBES_H

/*! \mainpage LibForBES documentation
 *
 * \section firstSteps Getting started
 *
 * ForBES is a solver for smooth and non-smooth convex optimization problems. 
 * LibForBES is a C++ implementation of ForBES powered the high-performance 
 * computational routines of BLAS, LAPACK and SuiteSparse.
 *
 * \section install_sec Installation
 *
 * LibForBES is available as a static library which you can easily use in your
 * C++ code.
 *  
 *
 */




/*! \page install Install and use libforbes
 * 
 * \brief First steps (install and start using)
 * 
 * \tableofcontents
 * 
 * Installing and using libforbes is quite easy.
 *
 * 
 */



/* Matrices */

/*! \page doc-matrix How to use Matrices
 * 
 * \brief Matrices and %Matrix operations in libforbes
 * 
 * \tableofcontents
 * Here we provide some examples on how to make the most out of the class
 * Matrix.
 * \section mat-sec0 Working with matrices
 * \subsection subsection1 First steps
 * Create an empty matrix and access its elements:
 * 
 * \code{.cpp} 
 * size_t n = 5, m = 7;
 * 
 * // Create an empty dense matrix 5-by-7:
 * Matrix A(n, m, Matrix::MATRIX_DENSE); 
 * 
 * // Set A(1,0) to a value:
 * A.set(1, 0, 10.56);
 * 
 * // retrieve the value of A(1,0):
 * double value = A.get(1, 0); 
 * \endcode
 * 
 * 
 * But you can also very easily create a random matrix:
 * 
 * \code{.cpp}
 * size_t n = 8;
 * 
 * // Make a random symmetric matrix:
 * Matrix S = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
 * 
 * // Make a random lower triangular matrix and random entries in [1, 3]:
 * Matrix L = MatrixFactory::MakeRandomMatrix(n, n, 1.0, 2.0, Matrix::MATRIX_LOWERTR);
 * 
 * // Make a random sparse matrix with given number of non-zero entries:
 * size_t nnz = 10;
 * Matrix Y = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.01, 1.5);
 * \endcode
 * 
 * You can very easily access the dimensions of a matrix:
 * \code{.cpp}
 * Matrix A;
 * size_t nrows = A.getNrows();
 * size_t ncols = A.getNcols();
 * \endcode
 *  
 * 
 * \subsection matops Matrix operations
 * Matrix addition, subtraction and multiplication is supported for all types of
 * matrices. You can check the type of the produced matrix using <code>getType</code>.
 * For instance, when adding or multiplying two sparse matrices, the result is a sparse matrix and
 * when adding two symmetric matrices, the result is a symmetric matrix.
 * Here is an example:
 * 
 * \code{.cpp}
 * size_t n = 30;
 * size_t nnz = 80;
 * Matrix A = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.01, 1.5);
 * Matrix B = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.01, 1.5);
 * Matrix C = A + B;
 * Matrix D = A - B;
 * \endcode
 * 
 * The operators <code>+=</code> and <code>-=</code> can also be used.
 * 
 * \code{.cpp}
 * // Matrix Multiplication
 * size_t n = 6;      // rows of A
 * size_t m = 4;      // cols of A, rows of B
 * size_t k = 5;      // cols of B
 * size_t nnz_A = 9;  // no. non-zero entries in A
 * size_t nnz_B = 7;  // no. non-zero entries in B
 * 
 * Matrix A = MatrixFactory::MakeRandomSparse(n, m, nnz_A, 0.01, 1.5);
 * Matrix B = MatrixFactory::MakeRandomSparse(m, k, nnz_B, 0.01, 1.5);
 * 
 * // Matrix C is n-by-k sparse:
 * Matrix C = A * B;
 * \endcode
 * 
 * You can also perform scalar multiplications
 * 
 * \code{.cpp}
 * Matrix A = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
 * A *= 3.465;              // A = A * 3.465
 * Matrix B = 7.0 * A;
 * Matrix C = 3.5 * B + A;
 * \endcode
 * 
 * 
 * 
 * \subsection qad Quadratic Forms
 * It may be often needed to compute quadratic forms
 * 
 * \f[
 *  F(x) = \frac{1}{2}x'Qx + q'x
 * \f]
 * 
 * To this end, two methods in <code>Matrix</code> have been implemented, both named
 * <code>Matrix::quad()</code>. Here is an example of use:
 * 
 * \code{.cpp}
 * const size_t n = 10;
 * Matrix Q = Matrix::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
 * Matrix x = Matrix::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
 * 
 * double f = Q.quad(x); // compute f = x'*Q*x
 * 
 * Matrix q = Matrix::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
 * 
 * double g = Q.quad(x, q); // compute f = (1/2)*x'*Q*x + q'*x
 * \endcode
 * 
 */

/* Matrix factorizations and solvers */

/*! \page doc-factorization Factorization of matrices
 *
 * \brief %Matrix factorizations and solvers
 * 
 * A matrix factorization is, in general, the procedure of decomposing a given 
 * matrix into a sum or product or other combination of two or more matrices which
 * possess favourable properties which can be exploited for instance to solve a 
 * linear system.
 * 
 * The ForBES API offers a simple interface to factor-solve procedures exporting 
 * to the client only what is necessary: a factor and a solve step. 
 * 
 * This is done by the abstract methods \c factorize() and \c solve() in FactoredSolver
 * which are implemented by derived classes.
 * 
 * \code
 * Matrix Q = ...;                      // A symmetric positive definite matrix
 * FactoredSolver * solver = new ...;   // Some factored solver
 * solver->factorize();                 // Factorize once
 * solver->solve(b, x);                 // Solve the system A*x = b
 * solver->solve(c, x);                 // Solve A*x = c
 * delete solver;                       // Release resources
 * \endcode
 * 
 * \tableofcontents 
 * 
 * \section Cholesky
 * The Cholesky factorization can be applied to any symmetric positive definite
 * matrix. 
 * 
 * A symmetric positive definite matrix \f$A\f$ can be written as \f$A=LL'\f$, where
 * \f$L\f$ is a lower-triangular matrix. 
 * 
 * Using this matrix one can solve linear systems of the form \f$Ax=b\f$ (for a given \f$b\f$)
 * at the cost of a forward and a backward substitution.
 * 
 * \section LDL-sec LDL factorization
 * The LDL factorization (often referred to as \f$LDL'\f$ can be applied to any
 * symmetric matrix. A given symmetric matrix \f$A\f$ is decomposed as \f$A=LDL'\f$
 * where \f$L\f$ is lower-triangular and \f$D\f$ is a diagonal matrix.
 * 
 * This is evidently similar to the Cholesky factorization described
 * above, but does not require that the factorized matrix be positive definite. 
 * LDL can be applied to any symmetric matrix.
 * 
 */

/* Linear operators */

/*! \page doc-linops Linear Operators
 *
 * \brief Working with linear operators
 * 
 * \tableofcontents 
 * 
 * \section linop-sec Linear Operator
 * 
 * A linear operator \f$T:\mathbb{R}^n \to \mathbb{R}^m\f$ may not always be 
 * available in the explicit form of a 
 * matrix, but in some implicit form to which we don't have direct access. 
 * 
 * Think of a linear operator <code>T</code> as a black box function for which we only
 * know that it is linear and we can evaluate its value \f$T(x)\f$ at a given
 * vector <code>x</code> and, sometimes, the value \f$T^*(x^*)\f$ of its adjoint
 * operator (which is again defined implicitly).
 * 
 * You can extends the class LinearOperator to define your own linear operators.
 * You will then need to implement the following methods:
 * 
 * \code{.cpp}
 * virtual Matrix call(Matrix& x) = 0;              // define T(x)
 * virtual Matrix callAdjoint(Matrix& x) = 0;       // define T^*(y)
 * virtual bool isSelfAdjoint() = 0;                // is it a self-adjoint operator?
 * virtual size_t dimensionIn() = 0;                // x - dimension
 * virtual size_t dimensionOut() = 0;               // T(x) - dimension
 * \endcode
 */

/* FUNCTIONS*/

/*! \page doc-functs Functions
 *
 * \brief ForBES functions
 * 
 * \tableofcontents 
 * 
 * \section functions-sec Functions
 * 
 * The ForBES function API blah blah
 * 
 * \subsection quad-fun-sec Quadratic
 * A quadratic function is a function of the form
 * 
 * \f[
 * F(x) = \frac{1}{2}x'Qx + q'x,
 * \f]
 * 
 * where <code>Q</code> is a symmetric positive definite matrix whose conjugate is given by
 * 
 * \f[
 * F^*(x^*) = (x^*-q)'Q^{-1}(x^*-q).
 * \f]
 * 
 * Here is a simple example:
 * \code
 * size_t n = 10;
 * Matrix Q = MatrixFactory::MakeRandomSparse(n, n, 20, 0.0, 1.0);
 * Matrix Eye = MatrixFactory::MakeIdentity(n, 10.0);
 * 
 * Matrix Qt(Q);        // Qt= Q
 * Qt.transpose();      // Qt = Qt'
 * Q += Qt;             // Q = Q + Qt (this will make Q a symmetric matrix)
 * 
 * Q += Eye;            // with this we ensure Q is positive definite
 *
 * Function *F = new Quadratic(Q); // F(x) = 0.5*x'Qx
 * Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
 * 
 * double f;            // value F(x)
 * double f_star;       // value F*(x) (conjugate of F at x)
 * Matrix grad;         // gradient of F* and x
 * 
 * int status = F->call(x, f);
 * status = F->callConj(x, f_star, grad);
 * 
 * std::cout << grad;   // print out the gradient
 * \endcode
 * 
 * The invocation of <code>callConj</code> involves the computation of a Cholesky
 * factor of <code>Q</code> which is stored internally in the instance of our 
 * quadratic function.
 * 
 * 
 * \subsection genquad-fun-sec Generalized Quadratic 
 * 
 * A generalized quadratic function is a function of the form
 * 
 * \f[
 * F(x) = \frac{1}{2}\langle x, T(x) \rangle
 * \f]
 * 
 * where <code>T</code> is a linear operator (here, an instance of LinearOperator).
 * This operator is expected to be self-adjoint.
 * 
 * 
 * Here is a very simple example of use where <code>T</code> is an instance of MatrixOperator.
 * 
 * \code
 * double fval;                         // the value of F at x, F(x) - to be computed
 * Matrix grad;                         // the gradient of F at x - to be computed
 * 
 * Matrix Q;
 * Matrix x;
 * 
 * LinearOperator *T = new MatrixOperator(Q);
 * Function *F = new QuadraticOperator(T);
 * 
 * int info = F -> call(x, fval, grad); // compute its value at x, F(x) and its gradient grad(F)(x)
 * 
 * delete T;
 * delete F;
 * \endcode
 * 
 * \subsection quadoveraff-fun-sec Quadratic-over affine
 * 
 * A quadratic-over-affine function is a function <code>F</code> is the form
 * 
 * \f[
 * F(x) = \frac{1}{2}x'Qx + q'x + \delta(x|E),
 * \f]
 * 
 * where <code>E</code> is an affine space
 * 
 * \f[
 * E = \{x: Ax = b\}
 * \f]
 * 
 * and delta is the indicator function
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
 * \f{align*}{
 * S g(x^*) &=  \begin{bmatrix}x^*-q\\b\end{bmatrix}\\
 * F^*(x^*) &= -\frac{1}{2} \left(g(x^*)'Qg(x^*) + (q-x^*)'g(x^*)\right)
 * \f}
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
 * 
 * \subsection indbox Indicator functions
 * There are four types of indicator functions implemented in ForBES, namely
 * (i) the indicator of a rectangle, that is a set of the form \f$B_{[l,u]}=\{x: l \leq x \leq u\}\f$,
 * (ii) the indicator of a set \f$B_{\geq l}=\{x: x \geq l\}\f$, (iii) or of a set \f$B_{\leq u}=\{x: x \leq u\}\f$
 * and (iv) the indicator of a ball \f$B_c=\{x: \|x\| \leq c\}\f$.
 * 
 * Here is an example:
 * 
 * \code{.cpp}
 * double lb = -1.0;    // uniform lower bound (-1.0)
 * double ub = 4.0;     // uniform upper bound (4.0)
 * 
 * Function *F = new IndBox(lb, ub);    // define the indicator function of a box
 * Matrix x(2, 1);
 * x[0] = -1.0;
 * x[1] = 4.0;
 * 
 * double fval;         // value F(x)
 * assert(ForBESUtils::STATUS_OK == F->call(x, fval));  // compute F(x)
 * assert(fval == 1.0);                                 // here F(x) = 1.0
 * 
 * x[0] = -1.1;
 * assert(ForBESUtils::STATUS_OK == F->call(x, fval));  // compute F(x)
 * assert(isinf(fval));                                 // here F(x) = infinity
 * 
 * delete F;
 * \endcode
 * 
 * One can use IndBox::callConj to compute the conjugate of the indicator function
 * of a box which is given by:
 * 
 * \f[
 * \delta^*(x^*|B_{[l,u]}) = \mathrm{mid}(x^*; l, u) = \min (\max(x^*, l), u)
 * \f]
 * 
 * 
 */

/* GROUPS */

/** \defgroup Functions ForBES Functions
 *
 * \brief ForBES functions
 * 
 * Lalala
 */


/** \defgroup Matrix Matrix
 *
 * \brief Matrices and utilities
 * 
 */


/*
 * This is a header you may include to your project - it will include all those
 * headers that are necessary to compile your project with libforbes. 
 * 
 */
#include "ForBESUtils.h"        /* ForBES utilities */

#include "Matrix.h"             /* Matrices */
#include "MatrixFactory.h"      /* Matrix Factory to construct matrices */

#include "LinearOperator.h"
#include "MatrixOperator.h"

#include "Function.h"           /* The Function API */
#include "Quadratic.h"          /* Quadratic functions */
#include "QuadOverAffine.h"
#include "QuadraticOperator.h"



#ifdef FORBES_TEST_UTILS        /* Define FORBES_TEST_UTILS in tests */

#ifndef TEST_UTILS_DEFINED
#define TEST_UTILS_DEFINED

#define _ASSERT_OK                          CPPUNIT_ASSERT_NO_THROW
#define _ASSERT                             CPPUNIT_ASSERT
#define _ASSERT_NOT(P)                      CPPUNIT_ASSERT(!(P))
#define _ASSERT_NUM_EQ(A,B,TOL)             CPPUNIT_ASSERT_DOUBLES_EQUAL((double)(A), (double)(B), (double)(TOL))
#define _ASSERT_EQ                          CPPUNIT_ASSERT_EQUAL
#define _ASSERT_NEQ(X,Y)                    CPPUNIT_ASSERT((X)!=(Y))
#define _ASSERT_EXCEPTION(P, EXCEPTION)     CPPUNIT_ASSERT_THROW(P, EXCEPTION)

#endif /* TEST_UTILS_DEFINED */
#endif /* FORBES_TEST_UTILS */

#endif	/* FORBES_H */

