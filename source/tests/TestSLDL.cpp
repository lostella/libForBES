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

    double beta = 1.13453;
    FactoredSolver * solver = new S_LDLFactorization(X, beta);
    int status = solver -> factorize();
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    status = solver -> solve(x, sol);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    Matrix g1 = X * Xt * sol;
    Matrix g2 = x - beta*sol;

    _ASSERT_EQ(g1, g2);

    delete solver;
}

void TestSLDL::testDenseShort() {
    const size_t n = 3;
    const size_t m = 7;
    Matrix A = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0); /* dense matrix*/
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    double beta = 0.5234234;
    Matrix sol;

    Matrix At(A);
    At.transpose();

    FactoredSolver * solver = new S_LDLFactorization(A, beta);
    int status = solver->factorize();
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    status = solver -> solve(x, sol);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);


    Matrix g1 = A * At * sol;
    Matrix g2 = x - beta*sol;

    delete solver;

    _ASSERT_EQ(g1, g2);

}

void TestSLDL::testDenseTall() {
    const size_t n = 10;
    const size_t m = 2;
    //Matrix X = MatrixFactory::MakeRandomSparse(n, m, 5, 0.0, 1.0);
    double X_data[20] = {
        0.004634224134067,
        0.774910464711502,
        0.817303220653433,
        0.868694705363510,
        0.084435845510910,
        0.399782649098896,
        0.259870402850654,
        0.800068480224308,
        0.431413827463545,
        0.910647594429523,
        0.181847028302852,
        0.263802916521990,
        0.145538980384717,
        0.136068558708664,
        0.869292207640089,
        0.579704587365570,
        0.549860201836332,
        0.144954798223727,
        0.853031117721894,
        0.622055131485066
    };
    Matrix X(n, m, X_data);

    double beta = 0.567;
    FactoredSolver * solver = new S_LDLFactorization(X, beta);
    int status = solver -> factorize();
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    Matrix x(10, 1);
    for (size_t i = 0; i < 10; i++) {
        x[i] = 1.0;
    }


    Matrix sol;
    status = solver -> solve(x, sol);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    double sol_expected_data[10] = {
        1.485739599580427,
        0.238313604427297,
        0.352734275537945,
        0.291783384355755,
        0.344089316703959,
        0.315246516232826,
        0.564132239648385,
        0.378781855470006,
        -0.138528061838890,
        -0.494174777885299,
    };
    Matrix sol_expected(10, 1, sol_expected_data);

    _ASSERT_EQ(sol_expected, sol);
    delete solver;

}

