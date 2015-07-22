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
                CPPUNIT_ASSERT(result.get(i, j) > 0);
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

void TestMatrixFactory::testMakeSparse() {
    int n = 5;
    int m = 8;
    int nnz = 10;
    Matrix *SP = new Matrix();
    for (int k = 0; k < 20; k++) {
        CPPUNIT_ASSERT_NO_THROW(*SP = MatrixFactory::MakeRandomSparse(n, m, nnz, 0.0, 1.0));
        CPPUNIT_ASSERT_EQUAL(n, SP->getNrows());
        CPPUNIT_ASSERT_EQUAL(m, SP->getNcols());
        CPPUNIT_ASSERT_EQUAL(Matrix::MATRIX_SPARSE, SP->getType());
        int actual_nnz = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (SP->get(i, j) != 0) {
                    actual_nnz++;
                }
            }
        }
        CPPUNIT_ASSERT_EQUAL(nnz, actual_nnz);
    }
    delete SP;
}

void TestMatrixFactory::testReadSparseFromFile() {
    FILE *fp;
    fp = fopen("sparse1.mx", "r");

    if (fp == NULL) {
        fprintf(stderr, "Can't open input file in.list!\n");
        exit(1);
    }

    Matrix A;
    CPPUNIT_ASSERT_NO_THROW(A = MatrixFactory::ReadSparse(fp));

    int n = 6;
    int m = 9;
    int nnz = 7;
    Matrix A_correct = MatrixFactory::MakeSparse(n, m, nnz, Matrix::SPARSE_UNSYMMETRIC);
    A_correct.set(0, 0, 3.8);
    A_correct.set(0, 1, 2.20);
    A_correct.set(1, 2, -1.18);
    A_correct.set(3, 3, 5.5);
    A_correct.set(3, 5, 1.23);
    A_correct.set(5, 1, 0.95);
    A_correct.set(5, 5, 2.68);
    
    CPPUNIT_ASSERT_EQUAL(A_correct, A);
}


void TestMatrixFactory::testSparse() {
    int n = 5;
    int m = 10;
    int max_nnz = 5;
    Matrix M = MatrixFactory::MakeSparse(n, m, max_nnz, Matrix::SPARSE_UNSYMMETRIC);
    CPPUNIT_ASSERT_EQUAL(Matrix::MATRIX_SPARSE, M.getType());
}

void TestMatrixFactory::testSparse2() {
    cholmod_common c;
    Matrix N = MatrixFactory::MakeSparse(20, 40, 6, Matrix::SPARSE_UNSYMMETRIC, &c);
    Matrix *M;
    CPPUNIT_ASSERT_NO_THROW(M = new Matrix(N));
    CPPUNIT_ASSERT_NO_THROW(delete M);
    CPPUNIT_ASSERT_EQUAL(0, c.status);
    CPPUNIT_ASSERT(c.memory_usage > 0);
    CPPUNIT_ASSERT_EQUAL(20, N.getNrows());
    CPPUNIT_ASSERT_EQUAL(40, N.getNcols());
}