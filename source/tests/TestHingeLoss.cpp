/*
 * File:   TestHingeLoss.cpp
 * Author: chung
 *
 * Created on Oct 29, 2015, 11:36:49 PM
 */

#include "TestHingeLoss.h"
#include "HingeLoss.h"
#include "MatrixFactory.h"
#include <iostream>


CPPUNIT_TEST_SUITE_REGISTRATION(TestHingeLoss);

TestHingeLoss::TestHingeLoss() {
}

TestHingeLoss::~TestHingeLoss() {
}

void TestHingeLoss::setUp() {
}

void TestHingeLoss::tearDown() {
}

void TestHingeLoss::testCall() {
    const size_t n = 140;
    Matrix b = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 10.0, Matrix::MATRIX_DENSE);
    const double mu = 1.5;
    Function * hinge = new HingeLoss(&b, mu);
    double f;
    int status = hinge->call(x, f);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    delete hinge;

}

void TestHingeLoss::testCallProx() {

    const size_t n = 5;
    const double bdata[] = {
        0.659605252908307,
        0.518594942510538,
        0.972974554763863,
        0.648991492712356,
        0.800330575352401
    };
    const double xdata[n] = {
        4.537977087269195,
        4.323915037834617,
        -8.253137954020456,
        0.834698148589140,
        1.331710076071617
    };
    Matrix b(n, 1, bdata);
    Matrix x(n, 1, xdata);

    const double mu = 0.7;
    const double gamma = 1.5;

    Function * hinge = new HingeLoss(&b, mu);
    Matrix prox(n, 1);
    hinge->callProx(x, gamma, prox);

    const double prox_expected_data[n] = {
        4.537977087269195,
        4.323915037834617,
        -6.791453098989512,
        0.648991492712356,
        1.331710076071617
    };
    Matrix prox_expected(n, 1, prox_expected_data);

    _ASSERT_EQ(prox_expected, prox);
    delete hinge;


}

