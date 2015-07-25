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
    double alpha = 2.5;
    Matrix result = MatrixFactory::MakeIdentity(n, alpha);
    _ASSERT_EQ(Matrix::MATRIX_DIAGONAL, result.getType());
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                _ASSERT_EQ(alpha, result.get(i, j));
            } else {
                _ASSERT_EQ(0.0, result.get(i, j));
            }
        }
    }

}

void TestMatrixFactory::testMakeRandomMatrix() {
    size_t n = 8;
    Matrix result = MatrixFactory::MakeRandomMatrix(n, n, 0.01, 1.5, Matrix::MATRIX_LOWERTR);
    _ASSERT_EQ(Matrix::MATRIX_LOWERTR, result.getType());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i < j) {
                _ASSERT_EQ(0.0, result.get(i, j));
            } else {
                _ASSERT(result.get(i, j) > 0);
            }
        }
    }

    result = MatrixFactory::MakeRandomMatrix(n, n, 0.01, 1.5, Matrix::MATRIX_SYMMETRIC);
    _ASSERT_EQ(Matrix::MATRIX_SYMMETRIC, result.getType());

    result = MatrixFactory::MakeRandomMatrix(n, n, 0.01, 1.5, Matrix::MATRIX_DIAGONAL);
    _ASSERT_EQ(Matrix::MATRIX_DIAGONAL, result.getType());

    result = MatrixFactory::MakeRandomMatrix(n, n, 0.01, 1.5, Matrix::MATRIX_DENSE);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, result.getType());

}

void TestMatrixFactory::testMakeSparse() {
    size_t n = 5;
    size_t m = 8;
    size_t nnz = 10;
    Matrix *SP = new Matrix();
    for (size_t k = 0; k < 20; k++) {
        _ASSERT_OK(*SP = MatrixFactory::MakeRandomSparse(n, m, nnz, 0.0, 1.0));
        _ASSERT_EQ(n, SP->getNrows());
        _ASSERT_EQ(m, SP->getNcols());
        _ASSERT_EQ(Matrix::MATRIX_SPARSE, SP->getType());
        size_t actual_nnz = 0;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < m; j++) {
                if (SP->get(i, j) != 0) {
                    actual_nnz++;
                }
            }
        }
        _ASSERT_EQ(nnz, actual_nnz);
    }
    delete SP;

    Matrix A;
    _ASSERT_OK(A = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0, Matrix::MATRIX_SPARSE));
    _ASSERT_EQ(Matrix::MATRIX_SPARSE, A.getType());
    _ASSERT_EQ(n, A.getNrows());
    _ASSERT_EQ(m, A.getNcols());
}

void TestMatrixFactory::testReadSparseFromFile() {
    FILE *fp;
    fp = fopen("matrices/sparse1.mx", "r");

    CPPUNIT_ASSERT_MESSAGE("File not found or cannot open", fp != NULL);

    Matrix A;
    _ASSERT_OK(A = MatrixFactory::ReadSparse(fp));

    _ASSERT_EQ(0, fclose(fp));

    size_t n = 6;
    size_t m = 9;
    size_t nnz = 7;
    Matrix A_correct = MatrixFactory::MakeSparse(n, m, nnz, Matrix::SPARSE_UNSYMMETRIC);
    A_correct.set(0, 0, 3.8);
    A_correct.set(0, 1, 2.20);
    A_correct.set(1, 2, -1.18);
    A_correct.set(3, 3, 5.5);
    A_correct.set(3, 5, 1.23);
    A_correct.set(5, 1, 0.95);
    A_correct.set(5, 5, 2.68);

    _ASSERT_EQ(A_correct, A);
}

void TestMatrixFactory::testSparse() {
    size_t n = 5;
    size_t m = 10;
    size_t max_nnz = 5;
    Matrix M = MatrixFactory::MakeSparse(n, m, max_nnz, Matrix::SPARSE_UNSYMMETRIC);
    _ASSERT_EQ(Matrix::MATRIX_SPARSE, M.getType());
}

void TestMatrixFactory::testSparse2() {
    const size_t n = 20;
    const size_t m = 40;
    const size_t nnz = 6;
    Matrix N = MatrixFactory::MakeSparse(n, m, nnz, Matrix::SPARSE_UNSYMMETRIC);
    Matrix *M;
    CPPUNIT_ASSERT_NO_THROW(M = new Matrix(N));
    CPPUNIT_ASSERT_NO_THROW(delete M);

    _ASSERT_EQ(n, N.getNrows());
    _ASSERT_EQ(m, N.getNcols());
    _ASSERT_EQ(0, Matrix::cholmod_handle()->status);
    _ASSERT(Matrix::cholmod_handle()->memory_usage > 0);
}