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

void TestOpDCT2::testAdj0() {
    const size_t n = 10;
    Matrix T(n, n);
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            T.set(k, i, std::cos(static_cast<double> (k) * M_PI * (static_cast<double> (i) + 0.5) / static_cast<double> (n)));
        }
    }

    T.transpose(); // defines the adjoint of DCT-II
    Matrix z = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0);

    Matrix Tz_correct = T*z;

    LinearOperator * dct2 = new OpDCT2(n);
    Matrix Tz = dct2->callAdjoint(z);
    _ASSERT_EQ(Tz_correct, Tz);

    delete dct2;
}

void TestOpDCT2::testCall1() {
    const size_t n = 10;
    Matrix T(n, n);
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            T.set(k, i, std::cos(static_cast<double> (k) * M_PI * (static_cast<double> (i) + 0.5) / static_cast<double> (n)));
        }
    }

    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0);

    Matrix y_correct = T*x;

    LinearOperator * dct2 = new OpDCT2(n);
    Matrix y = dct2->call(x);
    _ASSERT_EQ(y_correct, y);

    delete dct2;
}

void TestOpDCT2::testCall0() {

    const size_t n = 10;
    const double tol = 1e-6;

    LinearOperator * op = new OpDCT2(n);

    double x_data[n] = {
        0.368500000000000,
        -0.965700000000000,
        -1.332500000000000,
        1.484300000000000,
        -0.175200000000000,
        -0.373400000000000,
        -0.512200000000000,
        0.808700000000000,
        0.061900000000000,
        0.034600000000000
    };

    Matrix x(n, 1, x_data);

    double y_expected_data[n] = {
        -0.601000000000000,
        -1.162468863506907,
        -0.197505868217353,
        -0.411088627018780,
        0.384982166602636,
        4.028670175132235,
        2.343482143524826,
        -0.211797156859350,
        -0.624017833397364,
        -2.578437630903463
    };

    Matrix y_expected(n, 1, y_expected_data);

    Matrix y = op->call(x);
    _ASSERT_EQ(y_expected, y);


}

void TestOpDCT2::testCall() {
    const size_t n = 10;
    const size_t repeat = 50;
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
        Matrix Tstar_x = op->callAdjoint(*x);

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

        Matrix err = (*x) * T_y;
        Matrix temp = (*y) * Tstar_x;
        Matrix::add(err, -1.0, temp, 1.0);

        double error = std::abs(err.get(0, 0));
        _ASSERT(error < tol);
    }

    delete op;
    delete x;
    delete y;
}

void testOperatorLinearity(LinearOperator* op) {

    double a = 10.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
    double b = 10.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);

    Matrix x = MatrixFactory::MakeRandomMatrix(op->dimensionIn().first, op->dimensionIn().second, 0.0, 1.0);
    Matrix y = MatrixFactory::MakeRandomMatrix(op->dimensionIn().first, op->dimensionIn().second, 0.0, 1.0);

    // create z = ax + by
    Matrix z(x);
    Matrix::add(z, b, y, a);
    
    Matrix Tx = op->call(x);
    Matrix Ty = op->call(y);
    
    // create aTx + bTy
    Matrix T(Tx);
    Matrix::add(T, b, Ty, a);
    
    Matrix T2 = op->call(z);
    
    _ASSERT_EQ(T,T2);    
   
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




