/*
 * File:   TestIndBall2.cpp
 * Author: chung
 *
 * Created on Nov 3, 2015, 5:21:20 PM
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

#include "TestIndBall2.h"
#include "IndBall2.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestIndBall2);

TestIndBall2::TestIndBall2() {
}

TestIndBall2::~TestIndBall2() {
}

void TestIndBall2::setUp() {
}

void TestIndBall2::tearDown() {
}

void TestIndBall2::testCallProx() {
    double rho = 1.1;
    size_t n = 3;
    Matrix c(n, 1);
    c[0] = 1.0;
    c[1] = -0.5;
    c[2] = 3.0;

    Function * indB2 = new IndBall2(rho, c);

    Matrix x(n, 1);
    x[0] = 1.5;
    x[1] = 1.2;
    x[2] = 2.0;

    Matrix prox(n, 1);
    int status = indB2 -> callProx(x, 1.0, prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    const double tol = 1e-9;
    _ASSERT_NUM_EQ(1.270310252950645, prox[0], tol);
    _ASSERT_NUM_EQ(0.419054860032192, prox[1], tol);
    _ASSERT_NUM_EQ(2.459379494098711, prox[2], tol);

    double f_at_prox;
    status = indB2 -> callProx(x, 1.0, prox, f_at_prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_EQ(0.0, f_at_prox);

    delete indB2;
}

void TestIndBall2::testCategory() {
    Function * indB2 = new IndBall2();
    _ASSERT_NOT(indB2->category().defines_f());
    _ASSERT(indB2->category().defines_prox());
    _ASSERT_NOT(indB2->category().defines_conjugate());
    _ASSERT_NOT(indB2->category().defines_grad());
    delete indB2;
}

