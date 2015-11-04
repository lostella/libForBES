/*
 * File:   TestHingeLoss.cpp
 * Author: chung
 *
 * Created on Oct 29, 2015, 11:36:49 PM
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
    Function * hinge = new HingeLoss(b, mu);
    double f;
    _ASSERT(hinge->category().defines_f());
    int status = hinge->call(x, f);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    delete hinge;
}

void TestHingeLoss::testCall2() {
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

    Function * hinge = new HingeLoss(b, mu);

    double f;
    _ASSERT(hinge->category().defines_f());
    int status = hinge->call(x, f);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(6.641866860160479, f, 1e-11);
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

    Function * hinge = new HingeLoss(b, mu);
    Matrix prox(n, 1);
    _ASSERT(hinge->category().defines_prox());
    int status = hinge->callProx(x, gamma, prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    const double prox_expected_data[n] = {
        4.537977087269195,
        4.323915037834617,
        -6.791453098989512,
        0.648991492712356,
        1.331710076071617
    };
    Matrix prox_expected(n, 1, prox_expected_data);

    _ASSERT_EQ(prox_expected, prox);

    double f_at_prox;
    status = hinge->callProx(x, gamma, prox, f_at_prox);
    _ASSERT_EQ(prox_expected, prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(5.73070476830318, f_at_prox, 1e-12);

    delete hinge;
}

