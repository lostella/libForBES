/*
 * File:   TestNorm1.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Oct 30, 2015, 6:15:56 PM
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

#include "TestNorm1.h"
#include "Norm1.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestNorm1);

TestNorm1::TestNorm1() {
}

TestNorm1::~TestNorm1() {
}

void TestNorm1::setUp() {
}

void TestNorm1::tearDown() {
}

void TestNorm1::testCall() {
    const size_t n = 10;
    const double x_data[n] = {
        0.106216344928664,
        -0.372409740055537,
        0.198118402542975,
        -0.489687638016024,
        -0.339493413390758,
        0.951630464777727,
        0.920332039836564,
        -0.052676997680793,
        0.737858095516997,
        -0.269119426398556
    };
    const Matrix x(n, 1, x_data);
    double mu = 0.8;
    Function * norm1_fun = new Norm1(mu);
    double f;
    _ASSERT(norm1_fun->category().defines_f());
    int status = norm1_fun->call(const_cast<Matrix&>(x), f);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(3.550034050515676, f, 1e-12);
    delete norm1_fun;
}

void TestNorm1::testCallProx() {
    const size_t n = 10;
    const double x_data[n] = {
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
    Matrix x(n, 1, x_data);

    double mu = 0.8;
    Function * norm1_fun = new Norm1(mu);

    const double prox_expected_data[n] = {
        0,
        0,
        0,
        0.089687638016024,
        0,
        0.551630464777727,
        0.520332039836564,
        0,
        0.337858095516997,
        0
    };

    Matrix prox_expected(n, 1, prox_expected_data);

    const double gamma = 0.5;
    Matrix prox(n, 1);
    double value_at_prox;
    _ASSERT(norm1_fun->category().defines_prox());
    int status = norm1_fun->callProx(x, gamma, prox, value_at_prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    _ASSERT_EQ(prox_expected, prox);
    _ASSERT_NUM_EQ(1.199606590517850, value_at_prox, 1e-12);

    delete norm1_fun;
}

void TestNorm1::testCallProx2() {
    const size_t n = 10;
    const double x_data[n] = {
        0.537667139546100,
        1.833885014595086,
        -2.258846861003648,
        0.862173320368121,
        0.318765239858981,
        -1.307688296305273,
        -0.433592022305684,
        0.342624466538650,
        3.578396939725760,
        2.769437029884877
    };
    Matrix x(n, 1, x_data);

    double mu = 1.2;
    Function * norm1_fun = new Norm1(mu);

    const double prox_expected_data[n] = {
        0,
        1.233885014595086,
        -1.658846861003648,
        0.262173320368121,
        0,
        -0.707688296305273,
        0,
        0,
        2.978396939725760,
        2.169437029884877
    };

    Matrix prox_expected(n, 1, prox_expected_data);

    const double gamma = 0.5;
    Matrix prox(n, 1);
    _ASSERT(norm1_fun->category().defines_prox());
    int status = norm1_fun->callProx(x, gamma, prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    _ASSERT_EQ(prox_expected, prox);

    delete norm1_fun;
}

void TestNorm1::testDualNorm() {
    const double mu = 2.5;
    Norm1 normFun(mu);
    Matrix x(3, 1, Matrix::MATRIX_DENSE);
    x[0] = 0.5;
    x[1] = -0.9;
    x[2] = 0.85;

    double norm;
    double normDual;
    const double tol = 1e-12;

    int status = normFun.call(x, norm);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(5.6250, norm, tol);

    status = normFun.dualNorm(x, normDual);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(0.3600, normDual, tol);

    Norm1 normFun2;
    double normDual2;
    normFun2.dualNorm(x, normDual2);
    _ASSERT_NUM_EQ(normDual2 / mu, normDual, tol);
}

void TestNorm1::testConjugate() {
    const size_t n = 2;
    Matrix x(2,1); x[0] = 1;
    Norm1 norm;
    double f_1;
    double f_inf;
    double f_1_star;
    norm.call(x,f_1);
    norm.dualNorm(x, f_inf);
    _ASSERT_EQ(1.0, f_1);
    _ASSERT_EQ(1.0, f_inf);
    
    norm.callConj(x, f_1_star);
    _ASSERT_EQ(0.0, f_1_star);
    
    x[0] = 2.0;
    norm.callConj(x, f_1_star);
    _ASSERT(std::isinf(f_1_star));
    
    
}

