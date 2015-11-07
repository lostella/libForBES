/*
 * File:   TestOpDCT2.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Sep 15, 2015, 4:24:44 PM
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

#include "TestOpDCT2.h"
#include "OpAdjoint.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TestOpDCT2);

void testOperatorLinearity(LinearOperator* op);

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
    const double tol = 1e-7;

    LinearOperator * op = new OpDCT2(n);

    _ASSERT_EQ(n, op->dimensionIn().first);
    _ASSERT_EQ(n, op->dimensionOut().first);

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

void testOperatorLinearity(LinearOperator* op) {
    const size_t repeat = 300;
    const double tol = 1e-10;

    Matrix *x = new Matrix();
    Matrix *y = new Matrix();
    Matrix *ax = new Matrix();
    Matrix *by = new Matrix();
    Matrix *axby = new Matrix();
    Matrix *Taxby = new Matrix();
    Matrix *Tx = new Matrix();
    Matrix *Ty = new Matrix();
    Matrix *err = new Matrix();


    double a = 0.0;
    double b = 0.0;
    for (size_t r = 0; r < repeat; r++) {
        *x = MatrixFactory::MakeRandomMatrix(op->dimensionIn().first, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
        *y = MatrixFactory::MakeRandomMatrix(op->dimensionIn().first, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
        a = 10.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
        *ax = a * (*x);
        *by = b * (*y);
        *axby = *ax + *by;
        *Taxby = op->call(*axby);
        *Ty = op->call(*x);
        *Tx = op->call(*y);
        *err = *Taxby - a * (*Tx) - b * (*Ty);
        for (size_t j = 0; j < op->dimensionOut().first; j++) {
            _ASSERT(std::abs(err->get(j, 0)) < tol);
        }
    }

    delete x;
    delete y;
    delete ax;
    delete by;
    delete axby;
    delete Taxby;
    delete Tx;
    delete Ty;
    delete err;
}

void TestOpDCT2::testLinearity() {
    const size_t n = 50;
    LinearOperator * op = new OpDCT2(n);
    _ASSERT_EQ(n, op->dimensionIn().first);
    _ASSERT_EQ(n, op->dimensionOut().first);
    testOperatorLinearity(op);
    delete op;
}

void TestOpDCT2::testAdjointLinearity() {
    const size_t n = 100;
    LinearOperator * op = new OpDCT2(n);
    LinearOperator *adj = new OpAdjoint(*op);
    _ASSERT_EQ(n, adj->dimensionIn().first);
    _ASSERT_EQ(n, adj->dimensionOut().first);
    testOperatorLinearity(adj);
    delete op;
    delete adj;
}




