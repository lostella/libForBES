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

    _ASSERT_EQ(n, op->dimensionIn().first);
    _ASSERT_EQ(n - 1, op->dimensionOut().first);
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
    std::cout << "\nInvoking T* on y " << y.getNrows() << "x" << y.getNcols() << "\n";
    Matrix Tstar_y = op->callAdjoint(y);

    Matrix err = y * Tx;
    Matrix temp = x*Tstar_y;
    Matrix::add(err, -1.0, temp, 1.0);

    _ASSERT(std::abs(err.get(0, 0)) < tol);

    delete op;
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
