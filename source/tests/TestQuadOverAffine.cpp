/*
 * File:   TestQuadOverAffine.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Aug 27, 2015, 2:19:01 PM
 */

#include "TestQuadOverAffine.h"
#include "MatrixFactory.h"
#include "QuadOverAffine.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestQuadOverAffine);

TestQuadOverAffine::TestQuadOverAffine() {
}

TestQuadOverAffine::~TestQuadOverAffine() {
}

void TestQuadOverAffine::setUp() {
}

void TestQuadOverAffine::tearDown() {
}

void TestQuadOverAffine::testQuadOverAffine() {
    int n = 8;
    int s = 4;
    Matrix Q = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix A = MatrixFactory::MakeRandomMatrix(s, n, 0.0, -5.0, Matrix::MATRIX_DENSE);

    Matrix q = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix b = MatrixFactory::MakeRandomMatrix(s, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);

    QuadOverAffine * qoa;
    _ASSERT_OK(qoa = new QuadOverAffine(Q, q, A, b));
    _ASSERT_NEQ(NULL, qoa);

    Matrix y = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    double fstar = 0.0;
    Matrix grad;
    int status = qoa->callConj(y, fstar, grad);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NOT(std::abs(fstar) < 1e-7);

    _ASSERT_OK(delete qoa);

}

