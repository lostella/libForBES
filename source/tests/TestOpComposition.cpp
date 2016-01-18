/*
 * File:   TestOpComposition.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Sep 15, 2015, 2:40:14 AM
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

#include "TestOpComposition.h"
#include "OpReverseVector.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TestOpComposition);

TestOpComposition::TestOpComposition() {
}

TestOpComposition::~TestOpComposition() {
}

void TestOpComposition::setUp() {
}

void TestOpComposition::tearDown() {
}

void TestOpComposition::testCall() {
    const size_t n = 30;
    const double tol = 1e-8;
    Matrix A = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix B = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    MatrixOperator * T1 = new MatrixOperator(A); // T1(x) = A*x
    MatrixOperator * T2 = new MatrixOperator(B); // T2(x) = B*x

    OpComposition *G = new OpComposition(*T1, *T2); // G(x) = T1(T2(x))

    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE); // random vector x

    Matrix z = G -> call(x); // z = G(x) = T1(T2(x))
    Matrix t = T2 -> call(x); // t = T2(x)
    Matrix u = T1 -> call(t);
    Matrix err = z - u;
    for (size_t i = 0; i < n; i++) {
        _ASSERT_NUM_EQ(err.get(i, 0), 0.0, tol);
    }

    /* --- and now with different sizes */
    // x: m-by-1
    // A: p-by-n
    // B: n-ny-m
    const size_t m = 8;
    const size_t p = 15;
    A = MatrixFactory::MakeRandomMatrix(p, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    B = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0, Matrix::MATRIX_DENSE);

    T1 = new MatrixOperator(A);
    T2 = new MatrixOperator(B);


    G = new OpComposition(*T1, *T2);

    _ASSERT_EQ(p, T1->dimensionOut().first);
    _ASSERT_EQ(static_cast<size_t> (1), T1->dimensionOut().second);

    _ASSERT_EQ(p, G->dimensionOut().first);
    _ASSERT_EQ(static_cast<size_t> (1), G->dimensionOut().second);

    _ASSERT_EQ(m, G->dimensionIn().first);
    _ASSERT_EQ(static_cast<size_t> (1), G->dimensionIn().second);

    x = MatrixFactory::MakeRandomMatrix(m, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);

    z = G->call(x);
    t = T2 -> call(x); // t = T2(x)
    u = T1 -> call(t);
    err = z - u;
    for (size_t i = 0; i < G->dimensionOut().first; i++) {
        _ASSERT_NUM_EQ(err.get(i, 0), 0.0, tol);
    }

    delete T1;
    delete T2;
    delete G;

}

void TestOpComposition::testCall2() {
    size_t n = 10;
    size_t m = 15;

    LinearOperator * rev_op = new OpReverseVector(n);
    Matrix A = MatrixFactory::MakeRandomMatrix(m, n, -1.0, 2.0, Matrix::MATRIX_DENSE);
    LinearOperator * mat_op = new MatrixOperator(A);

    LinearOperator * op = new OpComposition(*mat_op, *rev_op);

    Matrix y = MatrixFactory::MakeRandomMatrix(m, 1, 0.0, 3.0, Matrix::MATRIX_DENSE);

    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 3.0, Matrix::MATRIX_DENSE);

    double alpha = -3.52540;
    double gamma = 1.41451;
    Matrix y_copy(y);
    int status = op->call(y, alpha, x, gamma);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    // verify result;
    Matrix x_rev = rev_op->call(x); // rev(x)
    Matrix t = A*x_rev; //
    t *= alpha;
    y_copy *= gamma;
    y_copy += t;
    _ASSERT_EQ(y_copy, y);

    delete rev_op;
    delete mat_op;
    delete op;
}

void TestOpComposition::testCallAdjoint() {
    const size_t p = 8;
    const size_t n = 15;
    const size_t m = 6;

    Matrix A = MatrixFactory::MakeRandomMatrix(p, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix B = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0, Matrix::MATRIX_DENSE);


    LinearOperator * A_op = new MatrixOperator(A);
    LinearOperator * B_op = new MatrixOperator(B);

    LinearOperator * AB_op = new OpComposition(*A_op, *B_op);

    Matrix x = MatrixFactory::MakeRandomMatrix(p, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);

    Matrix AB_op_adj_x = AB_op->callAdjoint(x);

    Matrix AB_op_adj_x_correct;
    B.transpose();
    A.transpose();
    AB_op_adj_x_correct = B * A * x;
    _ASSERT_EQ(AB_op_adj_x_correct, AB_op_adj_x);
    _ASSERT_NOT(AB_op->isSelfAdjoint());
    delete A_op;
    delete B_op;
    delete AB_op;

}

void TestOpComposition::testDimension() {
    const size_t n = 16;
    const size_t m = 9;
    const size_t p = 23;
    Matrix F = MatrixFactory::MakeRandomMatrix(p, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix A = MatrixFactory::MakeRandomMatrix(p, n + 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix B = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0, Matrix::MATRIX_DENSE);
    MatrixOperator *T1 = new MatrixOperator(A);
    MatrixOperator *T2 = new MatrixOperator(B);
    OpComposition *op;
    _ASSERT_EXCEPTION(op = new OpComposition(*T1, *T2), std::invalid_argument);

    T1 = new MatrixOperator(F);
    T2 = new MatrixOperator(B);

    op = new OpComposition(*T1, *T2);
    _ASSERT_EQ(m, op->dimensionIn().first);
    _ASSERT_EQ(p, op->dimensionOut().first);

    delete T1;
    delete T2;
    delete op;

}



