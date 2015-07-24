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

    Matrix M = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 10.0, Matrix::MATRIX_DENSE);

    LinearOperator *T = new MatrixOperator(M);

    Function *F = new QuadraticOperator(*T);
    Function *F2 = new Quadratic(M);

    double fval, fval2;
    const double tol = 1e-8;
    
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F2->call(x, fval2));
    _ASSERT_NUM_EQ(fval, fval2, tol);

    _ASSERT_OK(delete T);
    _ASSERT_OK(delete F);
    _ASSERT_OK(delete F2);
}

