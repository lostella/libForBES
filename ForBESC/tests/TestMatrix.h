/* 
 * File:   TestMatrix.h
 * Author: Chung
 *
 * Created on Jul 7, 2015, 8:07:04 PM
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


#ifndef TESTMATRIX_H
#define	TESTMATRIX_H

#include <cppunit/extensions/HelperMacros.h>

#define FORBES_TEST_UTILS
#include "ForBES.h"


class TestMatrix : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestMatrix);

    CPPUNIT_TEST(testMultiplication);
    CPPUNIT_TEST(testGetSet);
    CPPUNIT_TEST(testQuadratic);
    CPPUNIT_TEST(testQuadratic2);
    CPPUNIT_TEST(testAssignment);    
    CPPUNIT_TEST(testQuadratic3);    
    CPPUNIT_TEST(testAdditionBad);
    CPPUNIT_TEST(testAddition);
    CPPUNIT_TEST(testFBMatrix);
    CPPUNIT_TEST(testMakeRandomFBMatrix);
    CPPUNIT_TEST(testGetData);
    CPPUNIT_TEST(testGetNcols);
    CPPUNIT_TEST(testGetNrows);
    CPPUNIT_TEST(testIsColumnVector);
    CPPUNIT_TEST(testIsRowVector);
    CPPUNIT_TEST(testLength);
    CPPUNIT_TEST(testReshape);
    CPPUNIT_TEST(testReshapeBad);
    CPPUNIT_TEST(testQuadraticDot);
    CPPUNIT_TEST(testDiagonalGetSet);
    CPPUNIT_TEST(testDiagonalMultiplication);
    CPPUNIT_TEST(testDiagonalMultiplication2);
    CPPUNIT_TEST(testDenseTimesDiagonal);
    CPPUNIT_TEST(testQuadDiagonal);
    CPPUNIT_TEST(testQuadSymmetric);
    CPPUNIT_TEST(testSubtract);
    CPPUNIT_TEST(testCholesky);
    CPPUNIT_TEST(testSolveCholesky);
    CPPUNIT_TEST(testSolveCholeskyMulti);
    CPPUNIT_TEST(testLowerTriangular_getSet);
    CPPUNIT_TEST(testSymmetric_getSet);
    CPPUNIT_TEST(testSymmetricCholesky);
    CPPUNIT_TEST(testTranspose);
    CPPUNIT_TEST(testDiagonalTimesSymmetric);
    CPPUNIT_TEST(testDiagonalTimesLowerTri);
    CPPUNIT_TEST(testDenseTimesSymmetric);
    CPPUNIT_TEST(testDenseTimesLowerTriangular);
    CPPUNIT_TEST(testLeftTransposeMultiply);
    CPPUNIT_TEST(testRightTransposeMultiply);
    CPPUNIT_TEST(testLowerTriangularTraspose_getSet);
    CPPUNIT_TEST(testLeftSymmetricMultiply);
    CPPUNIT_TEST(testSparseGetSet);
    CPPUNIT_TEST(testSparseCholesky);
    CPPUNIT_TEST(testSparseDenseMultiply);
    CPPUNIT_TEST(testSparseSparseMultiply);
    CPPUNIT_TEST(testSparseAddDense);
    CPPUNIT_TEST(testSparseAddSparse);
    CPPUNIT_TEST(testSparseAddSparse2);
    CPPUNIT_TEST(testSparseQuad);
    CPPUNIT_TEST(testSparseQuadSparseX);
    CPPUNIT_TEST(testSparseQuad_q);
    CPPUNIT_TEST(testSparseDotProd);
    

    CPPUNIT_TEST_SUITE_END();

public:
    TestMatrix();
    virtual ~TestMatrix();
    void setUp();
    void tearDown();

private:
    void testMultiplication();
    void testQuadratic();
    void testQuadratic2();
    void testQuadratic3();
    void testQuadraticDot();
    void testGetSet();
    void testAssignment();
    void testAddition();
    void testAdditionBad();
    void testFBMatrix();
    void testMakeRandomFBMatrix();
    void testGetData();
    void testGetNcols();
    void testGetNrows();
    void testIsColumnVector();
    void testIsRowVector();
    void testLength();
    void testReshape();
    void testReshapeBad();    
    void testSubtract();    
    void testDiagonalGetSet();
    void testDiagonalMultiplication();
    void testDiagonalMultiplication2();
    void testDenseTimesDiagonal();
    void testQuadDiagonal();
    void testQuadSymmetric();
    void testCholesky();
    void testSolveCholesky();
    void testSolveCholeskyMulti();
    void testLowerTriangular_getSet();
    void testLowerTriangularTraspose_getSet();
    void testSymmetric_getSet();
    void testSymmetricCholesky();   
    void testDiagonalTimesSymmetric();
    void testDiagonalTimesLowerTri();
    void testDenseTimesSymmetric();
    void testDenseTimesLowerTriangular();    
    void testTranspose();
    void testLeftTransposeMultiply();
    void testRightTransposeMultiply();    
    void testLeftSymmetricMultiply();   
    void testSparseGetSet();
    void testSparseCholesky();
    void testSparseDenseMultiply();
    void testSparseSparseMultiply();
    void testSparseAddDense();
    void testSparseAddSparse();
    void testSparseAddSparse2();
    void testSparseQuad();
    void testSparseQuadSparseX();
    void testSparseQuad_q();
    void testSparseDotProd();
    
};

#endif	/* TESTMATRIX_H */

