/*
 * File:   TestOpComposition.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Sep 15, 2015, 2:40:14 AM
 */

#include "TestOpComposition.h"

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

    _ASSERT_EQ(p, T1->dimensionOut());
    _ASSERT_EQ(p, G->dimensionOut());
    _ASSERT_EQ(m, G->dimensionIn());

    x = MatrixFactory::MakeRandomMatrix(m, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);

    z = G->call(x);
    t = T2 -> call(x); // t = T2(x)
    u = T1 -> call(t);
    err = z - u;
    for (size_t i = 0; i < G->dimensionOut(); i++) {
        _ASSERT_NUM_EQ(err.get(i, 0), 0.0, tol);
    }

    delete T1;
    delete T2;
    delete G;

}

void TestOpComposition::testCallAdjoint() {
    const size_t p = 8;
    const size_t n = 15;
    const double tol = 1e-8;

    Matrix A = MatrixFactory::MakeRandomMatrix(p, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    LinearOperator * op = new MatrixOperator(A);

    Matrix y = MatrixFactory::MakeRandomMatrix(p, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);


    A.transpose();
    Matrix expected = A*y;
    A.transpose();

    Matrix op_star_x = op->callAdjoint(y);

    for (size_t i = 0; i < p; i++) {
        _ASSERT_NUM_EQ(expected.get(i, 0), op_star_x.get(i, 0), tol);
    }

    delete op;


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
    _ASSERT_EQ(m, op->dimensionIn());
    _ASSERT_EQ(p, op->dimensionOut());

    delete T1;
    delete T2;
    delete op;

}



