/* 
 * File:   ForBES.h
 * Author: Chung
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


/*
 * This is a header you may include to your project - it will include all those
 * headers that are necessary to compile your project with libforbes. 
 * 
 */
#include "ForBESUtils.h"        /* ForBES utilities */
#include "Matrix.h"             /* Matrices */
#include "MatrixFactory.h"      /* Matrix Factory to construct matrices */
#include "Function.h"           /* The Function API */
#include "Quadratic.h"          /* Quadratic functions */

#ifdef FORBES_TEST_UTILS        /* Define FORBES_TEST_UTILS in tests */

#ifndef TEST_UTILS_DEFINED
#define TEST_UTILS_DEFINED

#define _ASSERT_OK                          CPPUNIT_ASSERT_NO_THROW
#define _ASSERT                             CPPUNIT_ASSERT
#define _ASSERT_NOT(P)                      CPPUNIT_ASSERT(!(P))
#define _ASSERT_NUM_EQ(A,B,TOL)             CPPUNIT_ASSERT_DOUBLES_EQUAL((double)(A), (double)(B), (double)(TOL))
#define _ASSERT_EQ                          CPPUNIT_ASSERT_EQUAL
#define _ASSERT_NEQ(X,Y)                    CPPUNIT_ASSERT(X!=Y)
#define _ASSERT_EXCEPTION(P, EXCEPTION)     CPPUNIT_ASSERT_THROW(P, EXCEPTION)

#endif /* TEST_UTILS_DEFINED */
#endif /* FORBES_TEST_UTILS */

#endif	/* FORBES_H */

