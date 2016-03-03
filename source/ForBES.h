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

/*
 * This is a header you may include to your project - it will include all those
 * headers that are necessary to compile your project with libforbes. 
 * 
 */
#include "ForBESUtils.h"            /* ForBES utilities */

#include "FunctionOntologicalClass.h"

/*
 * MATRICES and FACTORIZATIONS
 */
#include "Matrix.h"                 /* Matrices */
#include "MatrixFactory.h"          /* Matrix Factory to construct matrices */
#include "LinSysSolver.h"           /* Abstraction tier for linear system solvers */
#include "FactoredSolver.h"         /* Generic factored solver tier */
#include "LDLFactorization.h"       /* LDL factorization */
#include "CholeskyFactorization.h"  /* Cholesky factorization */
#include "S_LDLFactorization.h"     /* LDL' factorization of AA'+bI */
#include "CGSolver.h"               /* Conjugate gradient solver (for linear operators) */
#include "MatrixSolver.h"           /* Factorized solver for matrices */

/* 
 * LINEAR OPERATORS
 */
#include "LinearOperator.h"         /* Generic linear operator (API) */
#include "MatrixOperator.h"         /* Matrix linear operator */
#include "OpAdjoint.h"              /* Adjoint of an operator */
#include "OpComposition.h"          /* Composition of linear operators */
#include "OpDCT2.h"                 /* Discrete Cosine Transform (DCT-II) */
#include "OpDCT3.h"                 /* Discrete Cosine Transform (DCT-III) */
#include "OpGradient.h"             /* Gradient of a vector and its conjugate */
#include "OpGradient2D.h"           /* 2D gradient (of matrices) */
#include "OpLTI.h"                  /* A linear time-invariant system */
#include "OpLinearCombination.h"    /* Linear combination of linear operators */
#include "OpReverseVector.h"        /* Vector reverse */
#include "OpSum.h"                  /* Sum of operators */


/*
 * FUNCTIONS
 */
#include "Function.h"                /* The Function API */
#include "Quadratic.h"               /* Quadratic functions */
#include "QuadOverAffine.h"          /* Quadratic over affine */
#include "QuadraticOperator.h"       /* Quadratic induced by a linear operator */
#include "DistanceToBox.h"           /* Distance-to-box */
#include "IndBox.h"                  /* Indicator of a box */ 
#include "IndPos.h"                  /* Indicator of a halfspace */
#include "IndSOC.h"                  /* Indicator of a second-order cone */
#include "IndBall2.h"                /* Indicator of a Euclidean ball */
#include "ElasticNet.h"              /* Elastic net regularization function */
#include "QuadraticLoss.h"           /* Quadratic loss function */
#include "QuadraticLossOverAffine.h" /* Quadratic loss over an affine subspace */
#include "LogLogisticLoss.h"         /* Log-logistic function */
#include "HingeLoss.h"               /* Hinge loss function */
#include "HuberLoss.h"               /* Huber loss function */
#include "Norm.h"                    /* Generic vector norm (abstract) */
#include "Norm1.h"                   /* Norm-1 */
#include "Norm2.h"                   /* Norm-2 */
#include "SumOfNorm2.h"              /* Sum of Norm-2 */
#include "SeparableSum.h"            /* Separable sum of proximable functions */
#include "ConjugateFunction.h"       /* Conjugate of a given function */

/*
 * FORBES SOLVER
 */
#include "IterativeSolver.h"         /* Generic iterative solver */
#include "FBCache.h"                 /* Low-level component of the libForBES engine */
#include "FBStopping.h"              /* Stopping criterion */
#include "FBStoppingRelative.h"      /* Stopping criterion using relative tolerance */
#include "FBProblem.h"               /* FB problem specifications */
#include "FBSplitting.h"             /* FB spliting algorithm */
#include "FBSplittingFast.h"         /* Accelerated FB splitting algorithm */


#ifdef FORBES_TEST_UTILS             /* Define FORBES_TEST_UTILS in tests */

#ifndef TEST_UTILS_DEFINED
#define TEST_UTILS_DEFINED

#define _ASSERT_OK                          CPPUNIT_ASSERT_NO_THROW
#define _ASSERT                             CPPUNIT_ASSERT
#define _ASSERT_NOT(P)                      CPPUNIT_ASSERT(!(P))
#define _ASSERT_NUM_EQ(A,B,TOL)             CPPUNIT_ASSERT_DOUBLES_EQUAL(static_cast<double>(A), static_cast<double>(B), static_cast<double>(TOL))
#define _ASSERT_EQ                          CPPUNIT_ASSERT_EQUAL
#define _ASSERT_NEQ(X,Y)                    CPPUNIT_ASSERT((X)!=(Y))
#define _ASSERT_EXCEPTION(P, EXCEPTION)     CPPUNIT_ASSERT_THROW(P, EXCEPTION)

#endif /* TEST_UTILS_DEFINED */
#endif /* FORBES_TEST_UTILS */

#endif	/* FORBES_H */

