/*
 * File:   TestOpDCT3.cpp
 * Author: chung
 *
 * Created on Sep 16, 2015, 3:42:02 AM
 */

#include "TestOpDCT3.h"
#include "OpAdjoint.h"
#include <cmath>

void testOperatorLinearity(LinearOperator*);

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
    _ASSERT_EQ(n, op->dimensionIn().first);
    _ASSERT_EQ(n, op->dimensionOut().first);
    Matrix *x = new Matrix();
    Matrix *y = new Matrix();

    for (size_t q = 0; q < repeat; q++) {
        *x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 10.0, Matrix::MATRIX_DENSE);
        *y = MatrixFactory::MakeRandomMatrix(n, 1, -3.0, 10.0, Matrix::MATRIX_DENSE);

        Matrix Tx = op->call(*x);
        Matrix Tstar_y = op->callAdjoint(*y);
        
        
        Matrix err = (*y) * Tx;
        Matrix temp = (*x) * Tstar_y;
        Matrix::add(err, -1.0, temp, 1.0);

        _ASSERT(std::abs(err.get(0, 0)) < tol);
    }

    delete op;
    delete x;
    delete y;
}

void testOperatorLinearity(LinearOperator* op) {

    double a = 10.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
    double b = 10.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);

    Matrix x = MatrixFactory::MakeRandomMatrix(op->dimensionIn().first, op->dimensionIn().second, 0.0, 1.0);
    Matrix y = MatrixFactory::MakeRandomMatrix(op->dimensionIn().first, op->dimensionIn().second, 0.0, 1.0);

    // create z = ax + by
    Matrix z(x);
    Matrix::add(z, b, y, a);
    
    Matrix Tx = op->call(x);
    Matrix Ty = op->call(y);
    
    // create aTx + bTy
    Matrix T(Tx);
    Matrix::add(T, b, Ty, a);
    
    Matrix T2 = op->call(z);
    
    _ASSERT_EQ(T,T2);    
   
}

void TestOpDCT3::testLinearity() {
    const size_t n = 50;
    LinearOperator * op = new OpDCT3(n);
    _ASSERT_EQ(n, op->dimensionIn().first);
    _ASSERT_EQ(n, op->dimensionOut().first);
    testOperatorLinearity(op);
    delete op;
}

void TestOpDCT3::testAdjointLinearity() {
    const size_t n = 60;
    LinearOperator * op = new OpDCT3(n);
    _ASSERT_EQ(n, op->dimensionIn().first);
    _ASSERT_EQ(n, op->dimensionOut().first);
    LinearOperator *adj = new OpAdjoint(*op);
    _ASSERT_EQ(n, adj->dimensionIn().first);
    _ASSERT_EQ(n, adj->dimensionOut().first);
    testOperatorLinearity(adj);
    delete op;
    delete adj;
}
