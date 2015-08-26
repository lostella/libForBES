/*
 * File:   TestLDL.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Aug 4, 2015, 8:18:33 PM
 */

#include "TestLDL.h"



CPPUNIT_TEST_SUITE_REGISTRATION(TestLDL);

TestLDL::TestLDL() {
}

TestLDL::~TestLDL() {
}

void TestLDL::setUp() {
}

void TestLDL::tearDown() {
}

void TestLDL::testSolveDense() {
    const size_t n = 40;
    const double tol = 1e-7;
    const size_t repetitions = 20;
    FactoredSolver *ldlSolver;

    for (size_t k = 0; k < repetitions; k++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
        Matrix At(A);
        At.transpose();
        A += At; // this will ensure that A is symmetric
        A *= 0.5;

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                _ASSERT_NUM_EQ(A.get(i, j), A.get(j, i), tol);
            }
        }

        ldlSolver = new LDLFactorization(A);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, ldlSolver->factorize());
        Matrix b = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
        Matrix x;
        _ASSERT_EQ(ForBESUtils::STATUS_OK, ldlSolver->solve(b, x));

        Matrix err = A * x - b;

        for (size_t i = 0; i < n; i++)
            _ASSERT(std::abs(err.get(i, 0)) < tol);

        delete ldlSolver;
    }
}

void TestLDL::testSolveSymmetric() {
    const size_t n = 8;
    const double tol = 1e-7;
    const size_t repetitions = 15;
    FactoredSolver *ldlSolver;
    Matrix b = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix x;
    for (size_t k = 0; k < repetitions; k++) {
        Matrix S = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
        ldlSolver = new LDLFactorization(S);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, ldlSolver->factorize());
        _ASSERT_EQ(ForBESUtils::STATUS_OK, ldlSolver->solve(b, x));
        Matrix err = S * x - b;
        for (size_t i = 0; i < n; i++)
            _ASSERT(std::abs(err.get(i, 0)) < tol);
        delete ldlSolver;
    }
}


