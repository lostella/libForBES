/*
 * File:   TestMatrixFactory.cpp
 * Author: Chung
 *
 * Created on Jul 13, 2015, 1:03:11 PM
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

#include "TestMatrixFactory.h"
#include "MatrixFactory.h"
#include <iostream>


CPPUNIT_TEST_SUITE_REGISTRATION(TestMatrixFactory);

TestMatrixFactory::TestMatrixFactory() {
}

TestMatrixFactory::~TestMatrixFactory() {
}

void TestMatrixFactory::setUp() {
}

void TestMatrixFactory::tearDown() {
}

void TestMatrixFactory::testMakeIdentity() {
    int n = 10;
    float alpha = 2.5;
    Matrix result = MatrixFactory::MakeIdentity(n, alpha);
    CPPUNIT_ASSERT_EQUAL(Matrix::MATRIX_DIAGONAL, result.getType());
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                CPPUNIT_ASSERT_EQUAL(alpha, result.get(i, j));
            } else {
                CPPUNIT_ASSERT_EQUAL(0.0f, result.get(i, j));
            }
        }
    }

}

void TestMatrixFactory::testMakeRandomMatrix() {
    int n = 8;
    Matrix result = MatrixFactory::MakeRandomMatrix(n, n, 0.01, 1.5, Matrix::MATRIX_LOWERTR);
    CPPUNIT_ASSERT_EQUAL(Matrix::MATRIX_LOWERTR, result.getType());
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i < j) {
                CPPUNIT_ASSERT_EQUAL(0.0f, result.get(i, j));
            } else {
                CPPUNIT_ASSERT(result.get(i, j)>0);
            }
        }
    }
    
    result = MatrixFactory::MakeRandomMatrix(n, n, 0.01, 1.5, Matrix::MATRIX_SYMMETRIC);
    CPPUNIT_ASSERT_EQUAL(Matrix::MATRIX_SYMMETRIC, result.getType());
    
    result = MatrixFactory::MakeRandomMatrix(n, n, 0.01, 1.5, Matrix::MATRIX_DIAGONAL);
    CPPUNIT_ASSERT_EQUAL(Matrix::MATRIX_DIAGONAL, result.getType());
    
    result = MatrixFactory::MakeRandomMatrix(n, n, 0.01, 1.5, Matrix::MATRIX_DENSE);
    CPPUNIT_ASSERT_EQUAL(Matrix::MATRIX_DENSE, result.getType());
    
}

