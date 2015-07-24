/*
 * File:   TestQuadraticOperator.cpp
 * Author: chung
 *
 * Created on Jul 24, 2015, 8:52:27 PM
 */

#include "TestQuadraticOperator.h"



CPPUNIT_TEST_SUITE_REGISTRATION(TestQuadraticOperator);

TestQuadraticOperator::TestQuadraticOperator() {
}

TestQuadraticOperator::~TestQuadraticOperator() {
}

void TestQuadraticOperator::setUp() {
}

void TestQuadraticOperator::tearDown() {
}

void TestQuadraticOperator::testCall() {
    size_t n = 10;

    Matrix Q = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix q = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 10.0, Matrix::MATRIX_DENSE);

    LinearOperator *T = new MatrixOperator(Q);
    LinearOperator *t = new MatrixOperator(q);

    Function *F = new QuadraticOperator(*T);
    Function *F2 = new Quadratic(Q);

    double fval, fval2;
    const double tol = 1e-8;
    
    Matrix grad, grad2;
    
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval, grad));
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F2->call(x, fval2, grad2));
    _ASSERT_NUM_EQ(fval, fval2, tol);
    _ASSERT_NOT(grad2.isEmpty());
    _ASSERT_NOT(grad.isEmpty());
    _ASSERT(grad.isColumnVector());
    _ASSERT(grad2.isColumnVector());
    _ASSERT_EQ(n, grad.getNrows());
    _ASSERT_EQ(grad, grad2);
       

    _ASSERT_OK(delete T);
    _ASSERT_OK(delete F);
    _ASSERT_OK(delete F2);
}

