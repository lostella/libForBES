/*
 * File:   TestOpReverseVector.cpp
 * Author: chung
 *
 * Created on Sep 16, 2015, 4:00:46 AM
 */

#include "TestOpReverseVector.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestOpReverseVector);

TestOpReverseVector::TestOpReverseVector() {
}

TestOpReverseVector::~TestOpReverseVector() {
}

void TestOpReverseVector::setUp() {
}

void TestOpReverseVector::tearDown() {
}

void TestOpReverseVector::testCall() {
    size_t n = 50;
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix y(x);
    LinearOperator *op = new OpReverseVector(n);
    Matrix xrev = op -> call(x);
    for (size_t i = 0; i < n; i++) {
        _ASSERT_EQ(y.get(i, 0), xrev.get(n - i - 1, 0));
    }
    delete op;
}

void TestOpReverseVector::testCallAdjoint() {
    size_t n = 80;
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix y = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    LinearOperator *op = new OpReverseVector(n);
    _ASSERT(op->isSelfAdjoint());

    Matrix Tx = op->call(x);
    Matrix Tstar_y = op->callAdjoint(y);

    Matrix err = y * Tx - x*Tstar_y;
    _ASSERT(std::abs(err.get(0, 0)) < 1e-10);
    delete op;
}


