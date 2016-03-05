/*
 * File:   TestLeastSquares.cpp
 * Author: chung
 *
 * Created on Mar 3, 2016, 2:47:16 PM
 */

#include "TestLeastSquares.h"
#include "MatrixFactory.h"
#include "LeastSquares.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestLeastSquares);

TestLeastSquares::TestLeastSquares() {
}

TestLeastSquares::~TestLeastSquares() {
}

void TestLeastSquares::setUp() {
}

void TestLeastSquares::tearDown() {
}

void TestLeastSquares::testSolveLS() {
    Matrix A = MatrixFactory::MakeRandomMatrix(5, 20, 0.0, 1.0);
    LeastSquares * ls = new LeastSquares(A);
    
    Matrix b = MatrixFactory::MakeRandomMatrix(5, 1, 1.0, 3.0);
    Matrix sol;
    ls->solve(b, sol);
    
    delete ls;
}


