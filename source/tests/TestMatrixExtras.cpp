/*
 * File:   TestMatrixExtras.cpp
 * Author: chung
 *
 * Created on Nov 8, 2015, 4:33:46 PM
 */

#include "TestMatrixExtras.h"
#include "MatrixFactory.h"
#include "ForBES.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestMatrixExtras);

TestMatrixExtras::TestMatrixExtras() {
}

TestMatrixExtras::~TestMatrixExtras() {
}

void TestMatrixExtras::setUp() {
}

void TestMatrixExtras::tearDown() {
}

void TestMatrixExtras::test_add_DD() {
    size_t n = 10;
    size_t m = 15;
    size_t repetitions = 200;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, m, 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(n, m, 2.0, 1.0);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        Matrix A_copy(A);
        Matrix B_copy(B);

        // A = gamma*A + alpha*B
        int status = Matrix::add(A, alpha, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
        _ASSERT_EQ(B_copy, B);

        A_copy *= gamma;
        B_copy *= alpha;
        Matrix C = A_copy + B_copy;

        _ASSERT_EQ(C, A);
    }

}

void TestMatrixExtras::test_add_SS() {
    size_t n = 10;
    size_t m = 15;
    size_t nnz = 100;
    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomSparse(n, m, nnz, 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomSparse(n, m, nnz, 2.0, 1.0);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        Matrix A_copy(A);
        Matrix B_copy(B);

        // A = gamma*A + alpha*B
        int status = Matrix::add(A, alpha, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
        _ASSERT_EQ(B_copy, B);

        A_copy *= gamma;
        B_copy *= alpha;
        Matrix C = A_copy + B_copy;

        _ASSERT_EQ(C, A);
    }
}


