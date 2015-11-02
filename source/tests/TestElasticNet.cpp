/*
 * File:   TestElasticNet.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Oct 29, 2015, 6:51:55 PM
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

#include "TestElasticNet.h"
#include "ElasticNet.h"
#include "ForBES.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestElasticNet);

TestElasticNet::TestElasticNet() {
}

TestElasticNet::~TestElasticNet() {
}

void TestElasticNet::setUp() {
}

void TestElasticNet::tearDown() {
}

void TestElasticNet::testCall() {
    Function * elastic = new ElasticNet(2.5, 1.3);
    const size_t n = 9;
    double xdata[n] = {-1.0, -3.0, 7.5, 2.0, -1.0, -1.0, 5.0, 2.0, -5.0};

    Matrix x(n, 1, xdata);

    double f;
    const double f_expected = 193.562500;
    const double tol = 1e-12;
    _ASSERT(elastic->category().defines_f());
    int status = elastic->call(x, f);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(f_expected, f, tol);
    delete elastic;
}

void TestElasticNet::testCallProx() {
    Function * elastic = new ElasticNet(2.5, 1.3);
    const size_t n = 9;
    double xdata[n] = {-1.0, -3.0, 7.5, 2.0, -1.0, -1.0, 5.0, 2.0, -5.0};
    const double gamma = 1.6;
    const double prox_expected_data[n] = {0.0, -0.1840, 1.0840, 0.0, 0.0, 0.0, 0.5840, 0.0, -0.5840};

    Matrix x(n, 1, xdata);
    Matrix prox_expected(n, 1, prox_expected_data);
    Matrix prox(n, 1);

    double f_at_prox;
    const double f_at_prox_expected = 5.5305800;
    const double tol = 1e-12;
    _ASSERT(elastic->category().defines_prox());
    int status = elastic->callProx(x, gamma, prox, f_at_prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(f_at_prox_expected, f_at_prox, tol);
    _ASSERT_EQ(prox_expected, prox);

    status = elastic->callProx(x, gamma, prox);
    _ASSERT_EQ(prox_expected, prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    
    delete elastic;
}

void TestElasticNet::testOther() {
    Function * elastic = new ElasticNet(2.5, 1.3);
    Matrix x;
    Matrix grad;
    double f_star;
    _ASSERT_NOT(elastic->category().defines_conjugate());
    int status = elastic->callConj(x, f_star);
    _ASSERT_EQ(ForBESUtils::STATUS_UNDEFINED_FUNCTION, status);

    _ASSERT_NOT(elastic->category().defines_conjugate_grad());
    status = elastic->callConj(x, f_star, grad);
    _ASSERT_EQ(ForBESUtils::STATUS_UNDEFINED_FUNCTION, status);
    delete elastic;
}
