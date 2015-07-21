/*
 * File:   TestQuadratic.cpp
 * Author: chung
 *
 * Created on Jul 9, 2015, 4:14:39 AM
 */

#include "TestQuadratic.h"
#include "Quadratic.h"

#include <cmath>

const static float MAT1[16] = {
    7, 2, -2, -1,
    2, 3, 0, -1,
    -2, 0, 3, -1,
    -1, -1, -1, 1
};
const static float MAT2[16] = {
    16.f, 2.f, 3.f, 13.f,
    5.f, 11.f, 10.f, 8.f,
    9.f, 7.f, 6.f, 12.f,
    4.f, 14.f, 15.f, 1.f
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestQuadratic);

TestQuadratic::TestQuadratic() {
}

TestQuadratic::~TestQuadratic() {
}

void TestQuadratic::setUp() {
}

void TestQuadratic::tearDown() {
}

void TestQuadratic::testQuadratic() {
    const float * Qdata;
    Qdata = MAT2;
    float xdata[4] = {-1.f, 1.f, 1.f, 1.f};

    Matrix Q = Matrix(4, 4, Qdata);
    Matrix x = Matrix(4, 1, xdata);
    Function *quad = new Quadratic(Q);

    float f;
    int info = quad -> call(x, f);
    CPPUNIT_ASSERT_EQUAL(Function::STATUS_OK, info);

    delete quad;
}

void TestQuadratic::testQuadratic2() {
    //  CPPUNIT_ASSERT(false);
}

void TestQuadratic::testQuadratic3() {
    //  CPPUNIT_ASSERT(false);
}

void TestQuadratic::testQuadratic4() {
    //   CPPUNIT_ASSERT(false);
}

void TestQuadratic::testCall() {
    const float * Qdata;
    Qdata = MAT2;
    float qdata[4] = {2.f, 3.f, 4.f, 5.f};
    float xdata[4] = {-1.f, 1.f, 1.f, 1.f};

    Matrix Q = Matrix(4, 4, Qdata);
    Matrix q = Matrix(4, 1, qdata);
    Matrix x = Matrix(4, 1, xdata);

    Quadratic quadratic(Q, q);
    float f = -999.0f;
    int status = quadratic.call(x, f);

    CPPUNIT_ASSERT_EQUAL(Function::STATUS_OK, status);
    CPPUNIT_ASSERT_EQUAL(74.0f, f);

    /* Second part */
    Quadratic quadratic2(Q);
    status = quadratic2.call(x, f);
    CPPUNIT_ASSERT_EQUAL(0, status);
    CPPUNIT_ASSERT_EQUAL(64.0f, f);
}

void TestQuadratic::testCallWithGradient() {    
    const int n = 4;
    const float * Qdata;
    Qdata = MAT1;
    float xdata[4] = {-1.f, 1.f, 1.f, 1.f};
    const float expected[4] = {-8, 0, 4, 0};

    Matrix Q = Matrix(4, 4, Qdata);
    Matrix x = Matrix(4, 1, xdata);

    Function * quad = new Quadratic(Q);
    float f = -999.0f;
    Matrix grad;
    CPPUNIT_ASSERT_EQUAL(Function::STATUS_OK, quad -> call(x, f, grad));

    for (int i = 0; i < n; i++) {
        CPPUNIT_ASSERT_EQUAL(expected[i], grad[i]);
    }

    delete quad;

}

void TestQuadratic::testCallConj() {
    const int n = 4;
    const float * Qdata;
    Qdata = MAT1;

    float qdata[4] = {2.f, 3.f, 4.f, 5.f};
    float xdata[4] = {-1.f, 1.f, 1.f, 1.f};

    Matrix Q = Matrix(4, 4, Qdata);
    Matrix q = Matrix(4, 1, qdata);
    Matrix x = Matrix(4, 1, xdata);

    Quadratic quadratic(Q, q);
    float fstar;
    quadratic.callConj(x, fstar);

    float expected = 421.0;
    CPPUNIT_ASSERT(std::fabs(expected - fstar) / expected < 1e-5);

    for (int i = 0; i < n; i++) {
        x[i] = 10 + 2 * i;
    }

    quadratic.callConj(x, fstar);
    expected = 3722;
    CPPUNIT_ASSERT(std::fabs(expected - fstar) / expected < 1e-5);
}

void TestQuadratic::testCategory() {
    Quadratic quadratic;
    int cat = quadratic.category();
    CPPUNIT_ASSERT_EQUAL(100, cat);
}

void TestQuadratic::testCallDiagonalMatrix() {
    const int n = 4;
    Matrix Q(n, n, Matrix::MATRIX_DIAGONAL);
    for (int i = 0; i < n; i++) {
        Q[i] = i + 2.0f;
    }

    Function *f = new Quadratic(Q);

    Matrix x(n, 1);
    for (int i = 0; i < n; i++) {
        x[i] = 2 * i + 1.0f;
    }

    float val;
    f -> call(x, val);

    CPPUNIT_ASSERT_EQUAL(374.f, val);
}
