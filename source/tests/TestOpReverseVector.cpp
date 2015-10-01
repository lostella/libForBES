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
    size_t n = 50;
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix y(x);
    LinearOperator *op = new OpReverseVector(n);
    Matrix xrev = op -> call(x);
    for (size_t i = 0; i < n; i++) {
        _ASSERT_EQ(y.get(i, 0), xrev.get(n - i - 1, 0));
    }
    delete op;
}

void TestOpReverseVector::testCallAdjoint() {
    size_t n = 80;
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix y = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    LinearOperator *op = new OpReverseVector(n);
    _ASSERT(op->isSelfAdjoint());

    Matrix Tx = op->call(x);
    Matrix Tstar_y = op->callAdjoint(y);

    Matrix err = y * Tx - x*Tstar_y;
    _ASSERT(std::abs(err.get(0, 0)) < 1e-10);
    delete op;
}


