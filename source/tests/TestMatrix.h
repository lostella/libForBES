/*
 * File:   TestMatrix.h
 * Author: Pantelis Sopasakis
 *
 * Created on Sep 14, 2015, 3:34:41 PM
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

#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <time.h>

class TestMatrix : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestMatrix);

    
    
    CPPUNIT_TEST(testGetSet);
    CPPUNIT_TEST(testGetSetTranspose);
    CPPUNIT_TEST(testQuadratic);
    CPPUNIT_TEST(testQuadratic2);
    CPPUNIT_TEST(testAssignment);
    CPPUNIT_TEST(testQuadratic3);
    CPPUNIT_TEST(testAdditionBad);
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
    CPPUNIT_TEST(testLowerTriangular_getSet);
    CPPUNIT_TEST(testSymmetric_getSet);
    CPPUNIT_TEST(testTranspose);
    CPPUNIT_TEST(testLeftTransposeMultiply);
    CPPUNIT_TEST(testRightTransposeMultiply);
    CPPUNIT_TEST(testLowerTriangularTraspose_getSet);
    CPPUNIT_TEST(testLeftSymmetricMultiply);
    CPPUNIT_TEST(testSparseGetSet);    
    CPPUNIT_TEST(testSparseAddDense);
    CPPUNIT_TEST(testSparseAddSparse);
    CPPUNIT_TEST(testSparseAddSparse2);
    CPPUNIT_TEST(testSparseQuad);
    CPPUNIT_TEST(testSparseQuadSparseX);
    CPPUNIT_TEST(testSparseQuad_q);
    CPPUNIT_TEST(testSparseDotProd);
    CPPUNIT_TEST(testSubmatrix);
    CPPUNIT_TEST(testSubmatrixSparse);
    CPPUNIT_TEST(testSubmatrixTranspose);
    CPPUNIT_TEST(testSubmatrixMultiply);
    CPPUNIT_TEST(testSubmatrixMultiplyTr);
    CPPUNIT_TEST(testToggleDiagonal);
    CPPUNIT_TEST(testOpplus);


    CPPUNIT_TEST(test_ADD1);
    CPPUNIT_TEST(test_ADD2);
    CPPUNIT_TEST(test_ADS);
    CPPUNIT_TEST(test_ADH);
    CPPUNIT_TEST(test_ADX);
    CPPUNIT_TEST(test_ADL);
    CPPUNIT_TEST(test_ADW);    

    CPPUNIT_TEST(test_ADDT);
    CPPUNIT_TEST(test_ADST);
    CPPUNIT_TEST(test_ADHT);
    CPPUNIT_TEST(test_ADXT);
    CPPUNIT_TEST(test_ADLT);
    CPPUNIT_TEST(test_ADWT);

    CPPUNIT_TEST(test_AHH);
    CPPUNIT_TEST(test_AHD);
    CPPUNIT_TEST(test_AHX);
    CPPUNIT_TEST(test_AHL);
    CPPUNIT_TEST(test_AHS);

    CPPUNIT_TEST(test_ASD);
    CPPUNIT_TEST(test_ASH);
    CPPUNIT_TEST(test_ASX);
    CPPUNIT_TEST(test_ASL);
    CPPUNIT_TEST(test_ASS);   
    
    CPPUNIT_TEST(test_ADTDT);
    CPPUNIT_TEST(test_ASTDT);
    
    CPPUNIT_TEST(test_ALL);
    CPPUNIT_TEST(test_ALX);

    CPPUNIT_TEST(test_AXX);
    
    CPPUNIT_TEST(test_ASST);
    
    

    CPPUNIT_TEST(test_EH);
    CPPUNIT_TEST(test_EX);
    CPPUNIT_TEST(test_EHT);
    CPPUNIT_TEST(test_EDT);
    CPPUNIT_TEST(test_EL);
    CPPUNIT_TEST(test_ES);
    CPPUNIT_TEST(test_EST);

    CPPUNIT_TEST(test_CD);
    CPPUNIT_TEST(test_CH);
    CPPUNIT_TEST(test_CS);

    CPPUNIT_TEST(test_MDD1);
    CPPUNIT_TEST(test_MDL);
    CPPUNIT_TEST(test_MDH);
    CPPUNIT_TEST(test_MDS);
    
    CPPUNIT_TEST(test_MSS);
    CPPUNIT_TEST(test_MSX);
    CPPUNIT_TEST(test_MSD);
    CPPUNIT_TEST(test_MSDT);
    CPPUNIT_TEST(test_MSTDT);
    
    CPPUNIT_TEST(test_MXH);
    CPPUNIT_TEST(test_MXL);
    CPPUNIT_TEST(test_MDX);

    CPPUNIT_TEST_SUITE_END();

        
    
public:
    TestMatrix();
    virtual ~TestMatrix();
    void setUp();
    void tearDown();

private:
    void testMethod();

    void testOpplus();
    void testQuadratic();
    void testQuadratic2();
    void testQuadratic3();
    void testQuadraticDot();
    void testGetSet();
    void testGetSetTranspose();
    void testAssignment();
    void test_ADD1();
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
    void testLowerTriangular_getSet();
    void testLowerTriangularTraspose_getSet();
    void testSymmetric_getSet();
    void testTranspose();
    void testLeftTransposeMultiply();
    void testRightTransposeMultiply();
    void testLeftSymmetricMultiply();
    void testSparseGetSet();    
    void test_MSS();
    void testSparseAddDense();
    void testSparseAddSparse();
    void testSparseAddSparse2();
    void testSparseQuad();
    void testSparseQuadSparseX();
    void testSparseQuad_q();
    void testSparseDotProd();
    void testSubmatrix();
    void testSubmatrixSparse();
    void testSubmatrixTranspose();
    void testSubmatrixMultiply();
    void testSubmatrixMultiplyTr();
    void testToggleDiagonal();

    /*
     * A: add
     * M: multiply
     * T: transpose
     * E: assignment
     * CL: scalar multiplication (left)
     * 
     * D: Dense
     * S: Sparse unsymmetric
     * H: Symmetric
     * X: Diagonal
     * L: Lower triangular
     * W: Sparse symmetric
     * 
     */

    /*
     * Copy-constructor and operator=
     */
    void test_EH();
    void test_EL();
    void test_EX();
    void test_ES();

    void test_EHT();
    void test_EDT();
    void test_EST();
    /*
     * Addition: X = A + B 
     */

    /*
     * DENSE + (?)
     */
    void test_ADD2();
    void test_ADS();
    void test_ADH();
    void test_ADX();
    void test_ADL();
    void test_ADW();


    /*
     * DENSE + (?)'
     */
    void test_ADDT();
    void test_ADST();
    void test_ADHT();
    void test_ADXT();
    void test_ADLT();
    void test_ADWT();
    
    /*
     * (?)' + DENSE'
     */
    void test_ADTDT();
    void test_ASTDT();

    /*
     * SYMMETRIC + (?) 
     */
    void test_AHD();
    void test_AHS();
    void test_AHH();
    void test_AHX();
    void test_AHL();
    void test_AHW();

    /*
     * SPARSE + (?)
     */
    void test_ASD();
    void test_ASS();
    void test_ASH();
    void test_ASX();
    void test_ASL();
    void test_ASW();


    /*
     * DIAGONAL + (?)
     */
    void test_AXD();
    void test_AXS();
    void test_AXH();
    void test_AXX();
    void test_AXL();
    void test_AXW();

    /*
     * SPARSE + (?)'
     */
    void test_ASST();

    /*
     * LOWER_TR + (?)
     */
    void test_ALX();
    void test_ALL();

    /*
     * alpha * X
     */
    void test_CD();
    void test_CH();
    void test_CL();
    void test_CS();
    void test_CX();

    /*
     * MULTIPLICATION TESTS
     */
    void test_MDD1();
    void test_MXH();
    void test_MXL();
    void test_MDH();
    void test_MDL();
    void test_MDX();
    void test_MSX();
    void test_MSD();
    void test_MDS();
    void test_MSDT();
    void test_MSTDT();
};

#endif	/* TESTMATRIX_H */

