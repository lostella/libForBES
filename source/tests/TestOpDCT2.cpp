/*
 * File:   TestOpDCT2.cpp
 * Author: chung
 *
 * Created on Sep 15, 2015, 4:24:44 PM
 */

#include "TestOpDCT2.h"
#include "LinearOperator.h"
#include "OpDCT2.h"
#include "MatrixFactory.h"


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
    LinearOperator * op = new OpDCT2();
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix y = op->call(x);    
    std::cout << y;
    delete op;
}



