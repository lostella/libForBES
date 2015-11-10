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

    Matrix A = MatrixFactory::MakeRandomSparse(n, m, nnz, 2.0, 1.0);
    for (size_t r = 0; r < repetitions; r++) {
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

void TestMatrixExtras::test_mult_DD() {
    size_t n = 8;
    size_t k = 6;
    size_t m = 5;
    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(k, m, 0.0, 1.0);
        Matrix C = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0);
        Matrix C_copy(C);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        int status = Matrix::mult(C, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        Matrix AB = A*B;
        AB *= alpha;
        C_copy *= gamma;
        C_copy += AB;

        _ASSERT_EQ(C_copy, C);
    }
}

void TestMatrixExtras::test_mult_DH() {
    size_t n = 8;
    size_t k = 6;
    size_t repetitions = 100;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(k, k, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
        Matrix C = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix C_copy(C);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        int status = Matrix::mult(C, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        Matrix AB = A*B;
        AB *= alpha;
        C_copy *= gamma;
        _ASSERT_NOT(C_copy == C);
        C_copy += AB;

        _ASSERT_EQ(C_copy, C);
    }
}

void TestMatrixExtras::test_mult_DDT() {
    size_t n = 8;
    size_t k = 6;
    size_t m = 5;
    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(m, k, 0.0, 1.0);
        Matrix C = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0);
        Matrix C_copy(C);

        B.transpose();

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        int status = Matrix::mult(C, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        Matrix AB = A*B;
        AB *= alpha;
        C_copy *= gamma;
        _ASSERT_NOT(C_copy == C);
        C_copy += AB;

        _ASSERT_EQ(C_copy, C);
    }
}

void TestMatrixExtras::test_mult_DX() {
    size_t n = 8;
    size_t k = 6;
    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(k, k, 0.0, 1.0, Matrix::MATRIX_DIAGONAL);
        Matrix C = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix C_copy(C);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        int status = Matrix::mult(C, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        Matrix AB = A*B;
        AB *= alpha;
        C_copy *= gamma;
        C_copy += AB;

        _ASSERT_EQ(C_copy, C);
    }
}

void TestMatrixExtras::test_mult_SS() {
    size_t n = 10;
    size_t k = 8;
    size_t m = 9;
    size_t nnz = 20;

    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {

        Matrix A = MatrixFactory::MakeRandomSparse(n, k, nnz, 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomSparse(k, m, nnz, 2.0, 1.0);
        Matrix C = MatrixFactory::MakeRandomSparse(n, m, nnz, 2.0, 1.0);

        double alpha = -2.0 + 2.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
        double gamma = -2.0 + 3.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);


        Matrix aAB = A*B;
        aAB *= alpha; /*     aAB = a * A * B          */

        Matrix gC(C);
        gC *= gamma; /*     gC  = g * C              */

        Matrix R = aAB + gC;

        Matrix C_copy(C);
        Matrix::mult(C_copy, alpha, A, B, gamma);

        _ASSERT_EQ(R, C_copy);
    }

}

void TestMatrixExtras::test_mult_SS2() {
    _ASSERT(false); // remember to build test with gamma = 0;
}


