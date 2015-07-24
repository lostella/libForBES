/*
 * File:   TestMatrixOperator.cpp
 * Author: chung
 *
 * Created on Jul 24, 2015, 8:39:54 PM
 */

#include "TestMatrixOperator.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestMatrixOperator);

TestMatrixOperator::TestMatrixOperator() {
}

TestMatrixOperator::~TestMatrixOperator() {
}

void TestMatrixOperator::setUp() {
    srand((unsigned)time(NULL));
}

void TestMatrixOperator::tearDown() {
}

void TestMatrixOperator::testCall() {
    size_t n = 10;
    size_t m = 3;
    
    Matrix M = MatrixFactory::MakeRandomMatrix(m, n, 0.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 10.0, Matrix::MATRIX_DENSE);
    LinearOperator *T = new MatrixOperator(M);
    Matrix y;
    _ASSERT_OK(y = T->call(x));
    
    Matrix z = M*x;
    _ASSERT_EQ(z, y);
    delete T;
    
}

