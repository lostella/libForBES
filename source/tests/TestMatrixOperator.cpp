/*
 * File:   TestMatrixOperator.cpp
 * Author: Pantelis Sopasakis
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
    srand(time(NULL));
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

void TestMatrixOperator::testCallAdjoint() {
    size_t n = 10;

    Matrix M = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 10.0, Matrix::MATRIX_SYMMETRIC);
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 10.0, Matrix::MATRIX_DENSE);
    MatrixOperator *T = new MatrixOperator(M);
    Matrix y;
    Matrix z;

    _ASSERT(T->isSelfAdjoint());
    _ASSERT_OK(y = T->call(x));
    _ASSERT_OK(z = T->callAdjoint(x));

    _ASSERT_EQ(z, y);


    M = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 10.0, Matrix::MATRIX_DENSE);
    _ASSERT_OK(T->SetMatrix(M));
    _ASSERT_NOT(T->isSelfAdjoint());
    
    M = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 10.0, Matrix::MATRIX_DENSE);
    _ASSERT_NOT(M.isSymmetric());
    _ASSERT_OK(T->SetMatrix(M));
    _ASSERT_NOT(T->isSelfAdjoint());
    _ASSERT_OK(y = T->call(x));
    _ASSERT_OK(z = T->callAdjoint(x));
    _ASSERT_NOT(y == z);
    
    delete T;
}

