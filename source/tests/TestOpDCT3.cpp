/*
 * File:   TestOpDCT3.cpp
 * Author: chung
 *
 * Created on Sep 16, 2015, 3:42:02 AM
 */

#include "TestOpDCT3.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TestOpDCT3);

TestOpDCT3::TestOpDCT3() {
}

TestOpDCT3::~TestOpDCT3() {
}

void TestOpDCT3::setUp() {
}

void TestOpDCT3::tearDown() {
}

void TestOpDCT3::testCall() {
    const size_t n = 15;
    const size_t repeat = 50;
    const double tol = 1e-10;

    LinearOperator * op = new OpDCT3(n);
    Matrix *x = new Matrix();
    Matrix *y = new Matrix();

    for (size_t q = 0; q < repeat; q++) {
        *x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 10.0, Matrix::MATRIX_DENSE);
        *y = MatrixFactory::MakeRandomMatrix(n, 1, -3.0, 10.0, Matrix::MATRIX_DENSE);

        Matrix Tx = op->call(*x);
        Matrix Tstar_y = op->callAdjoint(*y);
        Matrix err = (*y) * Tx - (*x) * Tstar_y;

        _ASSERT(std::abs(err.get(0, 0)) < tol);
    }

    delete op;
    delete x;
    delete y;
}

