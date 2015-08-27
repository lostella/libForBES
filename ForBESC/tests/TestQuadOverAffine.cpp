/*
 * File:   TestQuadOverAffine.cpp
 * Author: chung
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

void TestQuadOverAffine::testMethod() {
    int n = 8;
    int s = 4;
    Matrix Q = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix A = MatrixFactory::MakeRandomMatrix(s, n, 0.0, -5.0, Matrix::MATRIX_DENSE);

    Matrix q = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix b = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    
    QuadOverAffine * qoa = new QuadOverAffine(Q, q, A, b);

}

void TestQuadOverAffine::testFailedMethod() {
    //CPPUNIT_ASSERT(false);
}

