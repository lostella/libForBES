/*
 * File:   TestHuber.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Oct 30, 2015, 2:19:35 AM
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

#include "TestHuber.h"
#include "HuberLoss.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestHuber);

TestHuber::TestHuber() {
}

TestHuber::~TestHuber() {
}

void TestHuber::setUp() {
}

void TestHuber::tearDown() {
}

void TestHuber::testCall() {
    const size_t n = 10;
    const double xdata[n] = {
        0.106216344928664,
        0.372409740055537,
        0.198118402542975,
        0.489687638016024,
        0.339493413390758,
        0.951630464777727,
        0.920332039836564,
        0.052676997680793,
        0.737858095516997,
        0.269119426398556
    };
    Matrix x(n, 1, xdata);
    const double delta = 0.2;
    Function * huber = new HuberLoss(delta);

    double f;
    Matrix grad(n, 1);

    _ASSERT(huber->category().defines_f());
    int status = huber->call(x, f, grad);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    const double tol = 1e-12;
    _ASSERT_NUM_EQ(3.513800016594281, f, tol);
    
    status = huber->call(x, f);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(3.513800016594281, f, tol);

    const double grad_expected_data[n] = {
        0.531081724643320,
        1.000000000000000,
        0.990592012714875,
        1.000000000000000,
        1.000000000000000,
        1.000000000000000,
        1.000000000000000,
        0.263384988403965,
        1.000000000000000,
        1.000000000000000
    };

    Matrix grad_expected(n, 1, grad_expected_data);

    _ASSERT_EQ(grad_expected, grad);

    delete huber;
}

