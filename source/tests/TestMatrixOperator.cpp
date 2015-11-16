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

void TestMatrixOperator::testCall2() {
    size_t n = 10;
    size_t m = 3;
    Matrix A = MatrixFactory::MakeRandomMatrix(n, m, -2.0, 4.0);
    Matrix x = MatrixFactory::MakeRandomMatrix(m, 1, 0.0, 1.0);
    Matrix y = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0);

    Matrix y_copy(y);
    LinearOperator * op = new MatrixOperator(A);

    double alpha = M_PI;
    double gamma = M_SQRT2;

    int status = op->call(y, alpha, x, gamma); // y = gamma * y + alpha * A * x
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    y_copy *= gamma;
    Matrix z = A * x;
    z *= alpha;
    y_copy += z;

    _ASSERT_EQ(y_copy, y);

    delete op;
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

void TestMatrixOperator::testCallId() {
    size_t n = 10;
    Matrix Id = MatrixFactory::MakeIdentity(n, 1.0);
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1.0, 2.0, 2.0);
    MatrixOperator IdOp(Id);
    Matrix y = IdOp.call(x);
    _ASSERT_EQ(x, y);
}

void TestMatrixOperator::testCallAdjoint() {
    size_t n = 10;

    Matrix M = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 10.0, Matrix::MATRIX_SYMMETRIC);
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 10.0, Matrix::MATRIX_DENSE);
    MatrixOperator *T = new MatrixOperator(M);
    Matrix y;
    Matrix z;

    _ASSERT(T->isSelfAdjoint());
    y = T->call(x);
    z = T->callAdjoint(x);

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

