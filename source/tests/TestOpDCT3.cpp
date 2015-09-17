/*
 * File:   TestOpDCT3.cpp
 * Author: chung
 *
 * Created on Sep 16, 2015, 3:42:02 AM
 */

#include "TestOpDCT3.h"
#include "OpAdjoint.h"

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
    _ASSERT_EQ(n, op->dimensionIn());
    _ASSERT_EQ(n, op->dimensionOut());
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

void testOperatorLinearity(LinearOperator* op) {
    const size_t repeat = 300;
    const double tol = 1e-12;

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

void TestOpDCT3::testLinearity() {
    const size_t n = 50;
    LinearOperator * op = new OpDCT3(n);
    _ASSERT_EQ(n, op->dimensionIn());
    _ASSERT_EQ(n, op->dimensionOut());
    testOperatorLinearity(op);
    delete op;
}

void TestOpDCT3::testAdjointLinearity() {
    const size_t n = 60;
    LinearOperator * op = new OpDCT3(n);
    _ASSERT_EQ(n, op->dimensionIn());
    _ASSERT_EQ(n, op->dimensionOut());
    LinearOperator *adj = new OpAdjoint(*op);
    _ASSERT_EQ(n, adj->dimensionIn());
    _ASSERT_EQ(n, adj->dimensionOut());
    testOperatorLinearity(adj);
    delete op;
    delete adj;
}
