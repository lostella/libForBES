/*
 * File:   TestOpGradient.cpp
 * Author: chung
 *
 * Created on Sep 16, 2015, 1:57:49 AM
 */

#include "TestOpGradient.h"



CPPUNIT_TEST_SUITE_REGISTRATION(TestOpGradient);

TestOpGradient::TestOpGradient() {
}

TestOpGradient::~TestOpGradient() {
}

void TestOpGradient::setUp() {
}

void TestOpGradient::tearDown() {
}

void TestOpGradient::testCall() {
    const size_t n = 40;
    const double tol = 1e-8;

    LinearOperator * op = new OpGradient(n);

    _ASSERT_EQ(n, op->dimensionIn());
    _ASSERT_EQ(n - 1, op->dimensionOut());
    _ASSERT_NOT(op->isSelfAdjoint());
    
    Matrix x(n, 1);
    Matrix y(n - 1, 1);

    for (size_t i = 0; i < n; i++) {
        x.set(i, 0, i + 1);
    }

    for (size_t i = 0; i < n - 1; i++) {
        y.set(i, 0, 3 * i + 1);
    }

    Matrix Tx = op->call(x);
    Matrix Tstar_y = op->callAdjoint(y);

    Matrix err = y * Tx - x*Tstar_y;

    _ASSERT(std::abs(err.get(0, 0)) < tol);

    delete op;
}

