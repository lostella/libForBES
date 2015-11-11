/*
 * File:   TestOpReverseVector.cpp
 * Author: chung
 *
 * Created on Sep 16, 2015, 4:00:46 AM
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

#include "TestOpReverseVector.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestOpReverseVector);

TestOpReverseVector::TestOpReverseVector() {
}

TestOpReverseVector::~TestOpReverseVector() {
}

void TestOpReverseVector::setUp() {
}

void TestOpReverseVector::tearDown() {
}

void TestOpReverseVector::testCall() {
    for (size_t n = 5; n < 15; n++) {
        Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
        Matrix y(x);
        OpReverseVector op(n);
        _ASSERT_EQ(n, op.dimensionIn().first);
        _ASSERT_EQ(n, op.dimensionOut().first);
        Matrix xrev = op.call(x);
        for (size_t i = 0; i < n; i++) {
            _ASSERT_EQ(y.get(i, 0), xrev.get(n - i - 1, 0));
        }
    }

}

void TestOpReverseVector::testCallNotFixedSize() {
    size_t n1 = 50;
    size_t n2 = 50;
    Matrix x1 = MatrixFactory::MakeRandomMatrix(n1, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix x2 = MatrixFactory::MakeRandomMatrix(n2, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix y1(x1);
    Matrix y2(x2);
    LinearOperator *op = new OpReverseVector();
    size_t zero = 0;
    _ASSERT_EQ(zero, op->dimensionIn().first);
    _ASSERT_EQ(zero, op->dimensionOut().first);
    Matrix xrev1 = op -> call(x1);
    Matrix xrev2 = op -> call(x2);
    for (size_t i = 0; i < n1; i++) {
        _ASSERT_EQ(y1.get(i, 0), xrev1.get(n1 - i - 1, 0));
    }
    for (size_t i = 0; i < n2; i++) {
        _ASSERT_EQ(y2.get(i, 0), xrev2.get(n2 - i - 1, 0));
    }
    delete op;
}

void TestOpReverseVector::testCallAdjoint() {
    size_t n = 80;
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix y = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    LinearOperator *op = new OpReverseVector(n);
    _ASSERT(op->isSelfAdjoint());
    _ASSERT_EQ(n, op->dimensionIn().first);
    _ASSERT_EQ(n, op->dimensionOut().first);
    Matrix Tx = op->call(x);
    Matrix Tstar_y = op->callAdjoint(y);

    Matrix err = y * Tx;
    Matrix temp = x*Tstar_y;

    _ASSERT_EQ(err, temp);
    delete op;
}


