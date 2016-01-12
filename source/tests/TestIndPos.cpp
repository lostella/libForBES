/*
 * File:   TestIndPos.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Jan 12, 2016, 4:06:34 PM
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

#include <cmath>

#include "TestIndPos.h"
#include "IndPos.h"
#include "ForBES.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestIndPos);

TestIndPos::TestIndPos() {
}

TestIndPos::~TestIndPos() {
}

void TestIndPos::setUp() {
}

void TestIndPos::tearDown() {
}

void TestIndPos::testCall1() {
    Function * ind_pos = new IndPos();
    Matrix x(5, 1);
    x[1] = 10.0;
    double f = 123.0;
    int status = ind_pos -> call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(0.0, f, 1e-14);

    x[3] = -1;
    status = ind_pos -> call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));

    _ASSERT(std::isinf(f));

    delete ind_pos;
}

void TestIndPos::testCall2() {

    size_t n = 10;
    Matrix lb(n, 1);
    for (size_t i = 0; i < n; i++) {
        lb[i] = 3.0 + i;
    }

    Function * ind_pos = new IndPos(lb);

    Matrix x(lb);

    double f = 123.0;
    int status = ind_pos -> call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(0.0, f, 1e-14);

    x[3] = 4.9999;
    status = ind_pos -> call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT(std::isinf(f));

    delete ind_pos;
}

void TestIndPos::testCall3() {
    size_t n = 10;
    double lb = 1.50;

    Function * ind_pos = new IndPos(lb);

    Matrix x(n, 1);
    for (size_t i = 0; i < n; i++) {
        x[i] = lb;
    }

    double f = 123.0;
    int status = ind_pos -> call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(0.0, f, 1e-14);

    x[2] = 100.0;
    status = ind_pos -> call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(0.0, f, 1e-14);

    x[5] = 1.49999;
    status = ind_pos -> call(x, f);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT(std::isinf(f));

    delete ind_pos;
}

void TestIndPos::testConjugate1() {
    Function * ind_pos = new IndPos();
    size_t n = 200;
    Matrix y(n, 1);

    double f_star = 1234.567;
    int status = ind_pos->callConj(y, f_star);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(0.0, f_star, 1e-14);

    y[0] = -1.5;
    status = ind_pos->callConj(y, f_star);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(0.0, f_star, 1e-14);

    y[4] = 1.0;
    status = ind_pos->callConj(y, f_star);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT(std::isinf(f_star));

    delete ind_pos;
}

void TestIndPos::testConjugate2() {
    double uniform_lb = -1.50;
    Function * ind_pos = new IndPos(uniform_lb);
    size_t n = 200;
    Matrix y(n, 1); // y = 0

    double f_star;
    int status = ind_pos->callConj(y, f_star);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(0.0, f_star, 1e-14);



    y[0] = 1e-10;
    status = ind_pos->callConj(y, f_star);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT(std::isinf(f_star));

    y[0] = 1.0;
    y[1] = 2.0;
    status = ind_pos->callConj(y, f_star);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT(std::isinf(f_star));

    y[0] = -1.0;
    y[1] = -2.0;
    status = ind_pos->callConj(y, f_star);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(4.50, f_star, 1e-14);

    y[0] = -1.0;
    y[1] = -2.0;
    y[1] = 1.0;
    status = ind_pos->callConj(y, f_star);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT(std::isinf(f_star));

    delete ind_pos;
}

void TestIndPos::testConjugate3() {
    size_t n = 100;
    Matrix lb(n, 1);
    for (size_t i = 0; i < n; i++) {
        lb[i] = 3.0 + i;
    }
    Function * ind_pos = new IndPos(lb);

    Matrix y(n, 1);

    double f_star;
    int status = ind_pos->callConj(y, f_star);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(0.0, f_star, 1e-14);

    y[0] = 1.0;
    status = ind_pos->callConj(y, f_star);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT(std::isinf(f_star));

    y[0] = -1.0;
    y[1] = -2.0;
    y[2] = -3.0;
    status = ind_pos->callConj(y, f_star);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_NUM_EQ(-26.0, f_star, 1e-14);


    delete ind_pos;
}

void TestIndPos::testProx1() {
    Function * ind_pos = new IndPos();

    size_t n = 10;

    Matrix x(n, 1);
    x[0] = 10.0;
    x[1] = -10.0;
    x[n - 2] = -1.2;
    x[n - 1] = 3.4;

    Matrix z(n, 1);
    x.plusop(&z);

    Matrix prox(n, 1);
    int status = ind_pos -> callProx(x, 1.0, prox);

    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_EQ(z, prox);

    delete ind_pos;
}

void TestIndPos::testProx2() {
    size_t n = 50;
    Matrix lb(n, 1);
    Matrix prox(n, 1);
    Matrix x(n, 1);
    for (size_t i = 0; i < n; i++) {
        lb[i] = -3.0 + i;
        x[i] = lb[i] + std::pow(-1, i + 1);
    }
    Function * ind_pos = new IndPos(lb);

    int status = ind_pos -> callProx(x, 1.0, prox);
    _ASSERT(ForBESUtils::is_status_ok(status));

    for (size_t i = 0; i < n; i++) {
        if (i % 2 == 0) {
            _ASSERT_NUM_EQ(lb[i], prox[i], 1e-14);
        } else {
            _ASSERT_NUM_EQ(x[i], prox[i], 1e-14);
        }
    }

    delete ind_pos;
}

void TestIndPos::testProx3() {
    double uniform_lb = -1.50;
    Function * ind_pos = new IndPos(uniform_lb);
    size_t n = 10;
    Matrix x(n, 1); // y = 0
    

    Matrix prox(n, 1);    

    int status = ind_pos -> callProx(x, 1.0, prox);
    _ASSERT(ForBESUtils::is_status_ok(status));
    _ASSERT_EQ(x, prox);

    for (size_t i = 0; i < n; i++) {
        if (i % 2 == 0) {
            x[i] = -1.0*i - 2.0;
        } else {
            x[i] = 1.0*i + 2.0;
        }
    }
    status = ind_pos -> callProx(x, 1.0, prox);
    _ASSERT(ForBESUtils::is_status_ok(status));
    for (size_t i = 0; i < n; i++) {
        _ASSERT_NUM_EQ(std::max(uniform_lb, x[i]), prox[i], 1e-14);
        if (i % 2 == 0) {
            _ASSERT_NUM_EQ(uniform_lb, prox[i], 1e-14);
        } else {
            _ASSERT_NUM_EQ(i + 2.0, prox[i], 1e-14);
        }
    }

    delete ind_pos;
}

void TestIndPos::testCategory() {
    Function * ind_pos = new IndPos();
    FunctionOntologicalClass cat = ind_pos -> category();
    _ASSERT(cat.defines_conjugate());
    _ASSERT(cat.defines_prox());
    _ASSERT(cat.defines_f());
    _ASSERT_NOT(cat.defines_conjugate_grad());
    _ASSERT_NOT(cat.defines_grad());
    _ASSERT_NOT(cat.defines_hessian());
    _ASSERT_NOT(cat.defines_hessian_conj());
    delete ind_pos;
}


