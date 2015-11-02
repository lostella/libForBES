/*
 * File:   TestNorm2.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Oct 30, 2015, 9:53:14 PM
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

#include "TestNorm2.h"
#include "Norm2.h"
#include <iostream>


CPPUNIT_TEST_SUITE_REGISTRATION(TestNorm2);

TestNorm2::TestNorm2() {
}

TestNorm2::~TestNorm2() {
}

void TestNorm2::setUp() {
}

void TestNorm2::tearDown() {
}

void TestNorm2::testCall() {
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
    double mu = 1.1;
    double f;
    Function * norm2_fun = new Norm2(mu);
    
    _ASSERT(norm2_fun -> category().defines_f());
    int status = norm2_fun ->call(x, f);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(6.231255997317550, f, 1e-12);
    delete norm2_fun;
}

void TestNorm2::testCallProx() {
    const size_t n = 10;
    const double x_data[n] = {
        -5.399547760626085,
        12.139693865327418,
        2.901616899784422,
        -0.252219492758625,
        2.858971615304383,
        -0.819864233199099,
        -0.496577392865248,
        5.958790431141860,
        5.636137959201917,
        5.668769653718456
    };
    Matrix x(n, 1, x_data);
    Matrix prox(n, 1);
    double f_at_prox;
    double mu = 2.3760;
    Function * norm2_fun = new Norm2(mu);
    _ASSERT(norm2_fun -> category().defines_prox());
    double gamma = 0.50;
    _ASSERT(norm2_fun -> category().defines_f());
    int status = norm2_fun ->callProx(x, gamma, prox, f_at_prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(37.883327020592596, f_at_prox, 1e-12);

    const double prox_expected_data[n] = {
        -5.025125487612237,
        11.297887852628843,
        2.700409309223781,
        -0.234729769551485,
        2.660721188020410,
        -0.763012170143809,
        -0.462143094956132,
        5.545588445243951,
        5.245309749274328,
        5.275678655539919
    };
    Matrix prox_expected(n,1,prox_expected_data);
    _ASSERT_EQ(prox_expected, prox);
    
    status = norm2_fun ->callProx(x, gamma, prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_EQ(prox_expected, prox);
    
    x *= 0.05;
    status = norm2_fun ->callProx(x, gamma, prox, f_at_prox);
    _ASSERT_EQ(0.0, f_at_prox);
    
    delete norm2_fun;
}

void TestNorm2::testDualNorm() {
    const size_t n = 3;
    Matrix x(n,1);
    x[0]=1.0;
    x[1]=3.0;
    
    Norm2 norm;
    double f;
    double fd;
    norm.call(x,f);
    norm.dualNorm(x,fd);
    _ASSERT_NUM_EQ(3.162277660168380, f, 1e-12);
    _ASSERT_NUM_EQ(f, fd, 1e-10);
    
}
