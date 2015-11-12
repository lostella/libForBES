/*
 * File:   TestCholesky.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Aug 5, 2015, 12:28:55 AM
 * 
 * ForBES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *  
 * ForBES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with ForBES. If not, see <http://www.gnu.org/licenses/>.
 */

#include "TestCholesky.h"
#include <cmath>


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
    Matrix A_copy(A);
    Matrix b;

    CholeskyFactorization * cholFactorization = new CholeskyFactorization(A);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, cholFactorization ->factorize());

    // Goodbye A...
    A = Matrix();

    Matrix x;
    Matrix err;
    for (size_t k = 0; k < repetitions; k++) {
        b = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
        cholFactorization -> solve(b, x);
        err = A_copy * x - b;
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


    Matrix A_copy(A);
    Matrix x;
    FactoredSolver * solver = new CholeskyFactorization(A);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, solver -> factorize());
    A = Matrix(); // Goodbye A...
    _ASSERT_EQ(ForBESUtils::STATUS_OK, solver -> solve(b, x));

    Matrix ax = A_copy * x;
    for (size_t i = 0; i < n; i++) {
        _ASSERT(std::abs(ax.get(i, 0) - i - 1.0) < tol);
    }
    delete solver;
}

void TestCholesky::testCholeskySymmetric2() {
    const size_t n = 4;
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
    Matrix b(n, 1);
    for (size_t i = 0; i < n; i++) {
        b.set(i, 0, i + 1.0);
    }

    FactoredSolver * solver = new CholeskyFactorization(*A);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, solver -> factorize());
    
    // Bye bye Matrix A
    delete A;

    double sol_expected_data[4] = {-0.032619563872026,
        0.020584586456240,
        0.070701048816070,
        0.110212353853490};
    Matrix sol_expected(4, 1, sol_expected_data);

    Matrix sol;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, solver -> solve(b, sol));
    
    _ASSERT_EQ(sol_expected, sol);

    delete solver;
}

void TestCholesky::testCholeskySymmetric() {
    const double tol = 1e-7;
    size_t n = 30;
    size_t repetitions = 10;
    Matrix A = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
    for (size_t i = 0; i < n; i++) {
        A.set(i, i, A.get(i, i) + 1.2 * n);
    }
    Matrix A_copy(A);
    Matrix b;
    Matrix x;
    Matrix err;
    FactoredSolver * cholFactorization;
    cholFactorization = new CholeskyFactorization(A);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, cholFactorization ->factorize());

    A = Matrix();
    for (size_t k = 0; k < repetitions; k++) {
        b = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
        cholFactorization -> solve(b, x);
        err = A_copy * x - b;
        for (size_t i = 0; i < n; i++) {
            _ASSERT(std::abs(err.get(i, 0)) < tol);
        }
    }

    delete cholFactorization;
}



