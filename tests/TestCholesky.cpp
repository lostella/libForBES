/*
 * File:   TestCholesky.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Aug 5, 2015, 12:28:55 AM
 */

#include "TestCholesky.h"



CPPUNIT_TEST_SUITE_REGISTRATION(TestCholesky);

TestCholesky::TestCholesky() {
}

TestCholesky::~TestCholesky() {
}

void TestCholesky::setUp() {
}

void TestCholesky::tearDown() {
}

void TestCholesky::testCholeskyDense() {
    size_t n = 30;
    size_t repetitions = 80;
    const double tol = 1e-7;
    Matrix A = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix At = A;
    At.transpose();
    A += At;
    A *= 0.5;
    for (size_t i = 0; i < n; i++) {
        A.set(i, i, A.get(i, i) + 1.2 * n);
    }
    Matrix b;

    CholeskyFactorization * cholFactorization = new CholeskyFactorization(A);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, cholFactorization ->factorize());


    Matrix x;
    Matrix err;
    for (size_t k = 0; k < repetitions; k++) {
        b = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
        cholFactorization -> solve(b, x);
        err = A * x - b;
        for (size_t i = 0; i < n; i++) {
            _ASSERT(std::abs(err.get(i, 0)) < tol);
        }
    }
    delete cholFactorization;
}

void TestCholesky::testCholeskySparse() {
    const double tol = 1e-7;
    size_t n = 25;
    size_t nnz = 2 * n - 1;
    Matrix A = MatrixFactory::MakeSparseSymmetric(n, nnz);
    _ASSERT_EQ(Matrix::MATRIX_SPARSE, A.getType());
    Matrix b(n, 1);

    for (size_t i = 0; i < n; i++) {
        A.set(i, i, n + 2.5);
        b.set(i, 0, i + 1);
    }
    for (size_t i = 1; i < n; i++) { /* Set the LT part only */
        A.set(i, i - 1, 0.5);
    }


    Matrix x;
    FactoredSolver * solver = new CholeskyFactorization(A);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, solver -> factorize());
    _ASSERT_EQ(ForBESUtils::STATUS_OK, solver -> solve(b, x));

    Matrix ax = A * x;
    for (size_t i = 0; i < n; i++) {
        _ASSERT(std::abs(ax.get(i, 0) - i - 1.0) < tol);
    }
    delete solver;
}

void TestCholesky::testCholeskySymmetric2() {
    const int n = 4;
    Matrix *A = new Matrix(n, n, Matrix::MATRIX_SYMMETRIC);
    _ASSERT(A->isSymmetric());
    _ASSERT_EQ(Matrix::MATRIX_SYMMETRIC, A -> getType());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            A -> set(i, j, 3.2 * i + 0.2 * j + 0.45);
            if (i == j) A -> set(i, i, A->get(i, i) + 20.0);
        }
    }

    Matrix L;
    Matrix x(n, 1);
    for (size_t i = 0; i < n; i++) {
        x.set(i, 0, i + 1.0);
    }

    FactoredSolver * solver = new CholeskyFactorization(*A);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, solver -> factorize());


    delete A;
}

void TestCholesky::testCholeskySymmetric() {
    const double tol = 1e-7;
    size_t n = 30;
    size_t repetitions = 10;
    Matrix A = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
    for (size_t i = 0; i < n; i++) {
        A.set(i, i, A.get(i, i) + 1.2 * n);
    }
    Matrix b;
    Matrix x;
    Matrix err;
    FactoredSolver * cholFactorization;
    cholFactorization = new CholeskyFactorization(A);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, cholFactorization ->factorize());

    for (size_t k = 0; k < repetitions; k++) {
        b = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
        cholFactorization -> solve(b, x);
        err = A * x - b;
        for (size_t i = 0; i < n; i++) {
            _ASSERT(std::abs(err.get(i, 0)) < tol);
        }
    }

    delete cholFactorization;
}




