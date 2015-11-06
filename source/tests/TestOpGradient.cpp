/*
 * File:   TestOpGradient.cpp
 * Author: chung
 *
 * Created on Sep 16, 2015, 1:57:49 AM
 */

#include "TestOpGradient.h"
#include "OpAdjoint.h"


void testOperatorLinearity(LinearOperator* op);

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

void testOperatorLinearity(LinearOperator* op) {
    const size_t repeat = 50;
    const double tol = 1e-7;

    Matrix *x = new Matrix();
    Matrix *y = new Matrix();
    Matrix *ax = new Matrix();
    Matrix *by = new Matrix();
    Matrix *axby = new Matrix();
    Matrix *Taxby = new Matrix();
    Matrix *Tx = new Matrix();
    Matrix *Ty = new Matrix();
    Matrix *err = new Matrix();


    double a = 0.0;
    double b = 0.0;
    for (size_t r = 0; r < repeat; r++) {
        *x = MatrixFactory::MakeRandomMatrix(op->dimensionIn(), 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
        *y = MatrixFactory::MakeRandomMatrix(op->dimensionIn(), 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
        a = 10.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
        *ax = a * (*x);
        *by = b * (*y);
        *axby = *ax + *by;
        *Taxby = op->call(*axby);
        *Ty = op->call(*x);
        *Tx = op->call(*y);
        *err = *Taxby - a * (*Tx) - b * (*Ty);
        for (size_t j = 0; j < op->dimensionOut(); j++) {
            _ASSERT(std::abs(err->get(j, 0)) < tol);
        }
    }

    delete x;
    delete y;
    delete ax;
    delete by;
    delete axby;
    delete Taxby;
    delete Tx;
    delete Ty;
    delete err;
}

void TestOpGradient::testLinearity() {
    const size_t n = 50;
    LinearOperator * op = new OpGradient(n);
    testOperatorLinearity(op);
    delete op;
}

void TestOpGradient::testAdjointLinearity() {
    const size_t n = 50;
    LinearOperator * op = new OpGradient(n);
    LinearOperator *adj = new OpAdjoint(*op);
    testOperatorLinearity(adj);
    delete op;
    delete adj;
}
