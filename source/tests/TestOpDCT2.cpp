/*
 * File:   TestOpDCT2.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Sep 15, 2015, 4:24:44 PM
 */

#include "TestOpDCT2.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TestOpDCT2);

TestOpDCT2::TestOpDCT2() {
}

TestOpDCT2::~TestOpDCT2() {
}

void TestOpDCT2::setUp() {
}

void TestOpDCT2::tearDown() {
}

void TestOpDCT2::testCall() {
    const size_t n = 10;
    const size_t repeat = 100;
    const double tol = 1e-8;

    LinearOperator * op = new OpDCT2(n);

    _ASSERT_EQ(n, op->dimensionIn());
    _ASSERT_EQ(n, op->dimensionOut());

    _ASSERT_NOT(op->isSelfAdjoint());

    Matrix *x = new Matrix();
    Matrix *y = new Matrix();

    for (size_t i = 0; i < repeat; i++) {
        *x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
        *y = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);

        Matrix T_y = op->call(*y);
        Matrix Tstar_x = op->callAdjoint(*y);

        // make sure T_y and Tstar_x are not all zero
        bool T_y_OK = false;
        bool Tstar_x_OK = false;
        for (size_t j = 0; j < n; j++) {
            if (std::abs(T_y.get(j, 0)) > tol) {
                T_y_OK = true;
                break;
            }
        }
        for (size_t j = 0; j < n; j++) {
            if (std::abs(Tstar_x.get(j, 0)) > tol) {
                Tstar_x_OK = true;
                break;
            }
        }
        
        _ASSERT(T_y_OK);
        _ASSERT(Tstar_x_OK);
        
        Matrix err = (*x) * T_y - (*y) * Tstar_x;

        double error = std::abs(err.get(0, 0));
        _ASSERT(error < tol);
    }

    delete op;
    delete x;
    delete y;
}



