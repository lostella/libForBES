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

void TestLDL::testSolveSparse() {
    const double tol = 1e-7;
    const size_t N = 10;
    const size_t NNZ = 19;

    Matrix A = MatrixFactory::MakeSparse(N, N, NNZ, Matrix::SPARSE_SYMMETRIC_L);
    A.set(0, 0, 1.7);
    A.set(0, 8, 0.13);
    A.set(1, 1, 1.0);
    A.set(1, 9, 0.01);
    A.set(2, 2, 1.5);
    A.set(3, 3, 1.1);
    A.set(4, 4, 2.6);
    A.set(5, 5, 1.2);
    A.set(6, 6, 1.3);
    A.set(6, 6, 1.3);
    A.set(7, 7, 1.6);
    A.set(8, 8, 1.4);
    A.set(9, 9, 3.1);
    A.set(1, 4, 0.02);
    A.set(4, 6, 0.16);
    A.set(4, 7, 0.09);
    A.set(4, 8, 0.52);
    A.set(4, 9, 0.53);
    A.set(6, 9, 0.56);
    A.set(7, 8, 0.11);


    LDLFactorization * solver = new LDLFactorization(A);
    int status = solver->factorize();
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    double b[N] = {.287, .22, .45, .44, 2.486, .72, 1.55, 1.424, 1.621, 3.759};
    Matrix rhs(N, 1, b);

    Matrix sol;
    status = solver->solve(rhs, sol);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);


    for (size_t i = 0; i < N; i++)
        _ASSERT_NUM_EQ((i + 1.0) / 10.0, sol.get(i, 0), tol);

    Matrix err = A * sol - rhs;
    for (size_t i = 0; i < N; i++)
        _ASSERT(std::abs(err.get(i, 0)) < tol);

    delete solver;
}



