/*
 * File:   TestSLDL.cpp
 * Author: chung
 *
 * Created on Nov 5, 2015, 4:08:41 AM
 */

#include "TestSLDL.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestSLDL);

TestSLDL::TestSLDL() {
}

TestSLDL::~TestSLDL() {
}

void TestSLDL::setUp() {
}

void TestSLDL::tearDown() {
}

void TestSLDL::testFactorizeAndSolve() {
    const size_t n = 5;
    const size_t m = 2;
    Matrix X = MatrixFactory::MakeRandomSparse(n, m, 5, 0.0, 1.0);

    Matrix Xt(X);
    Xt.transpose();

    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix sol;

    double beta = 1;
    FactoredSolver * solver = new S_LDLFactorization(X, beta);
    int status = solver -> factorize();
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    status = solver -> solve(x, sol);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    
    Matrix g1 = X * Xt * sol;
    Matrix g2 = x - sol;

    _ASSERT_EQ(g1, g2);

    delete solver;
}

