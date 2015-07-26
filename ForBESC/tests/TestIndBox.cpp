/*
 * File:   TestIndBox.cpp
 * Author: chung
 *
 * Created on Jul 26, 2015, 5:41:14 PM
 */

#include "TestIndBox.h"



CPPUNIT_TEST_SUITE_REGISTRATION(TestIndBox);

TestIndBox::TestIndBox() {
}

TestIndBox::~TestIndBox() {
}

void TestIndBox::setUp() {
}

void TestIndBox::tearDown() {
}

void TestIndBox::testCall() {
    double lb = -1.0;
    double ub = 4.0;
    Function *F = new IndBox(lb, ub);

    int funCat = F->category();
    _ASSERT_EQ(Function::CAT_INDICATOR, funCat);

    Matrix x(2, 1);
    x[0] = -1.0;
    x[1] = 4.0;

    double fval = -1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_EQ(1.0, fval);

    x[1] = 5.0;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT(isinf(fval));

    x[0] = -0.5;
    x[1] = 3.1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_NOT(isinf(fval));
    _ASSERT_EQ(1.0, fval);

    x[0] = -1 - 1e-9;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_NEQ(1.0, fval);
    _ASSERT(isinf(fval));


    _ASSERT_OK(delete F);

}

void TestIndBox::testCall2() {
    const size_t n = 10;
    Matrix lb(n, 1);
    Matrix ub(n, 1);
    for (int i = 0; i < n; i++) {
        lb[i] = -i - 1;
        ub[i] = i + 1;
    }

    Function *F;
    _ASSERT_OK(F = new IndBox(lb, ub));

    Matrix x(n, 1);

    double fval = -1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_EQ(1.0, fval);

    x = ub;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_EQ(1.0, fval);

    x = lb;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_EQ(1.0, fval);


    Matrix y(n, n);
    _ASSERT_EXCEPTION(F->call(y, fval), std::invalid_argument);

    lb = Matrix(n, n);
    _ASSERT_EXCEPTION(F = new IndBox(lb, ub), std::invalid_argument);

    lb = Matrix(n, 1);
    ub = Matrix(n, n);

    _ASSERT_EXCEPTION(F = new IndBox(lb, ub), std::invalid_argument);
    _ASSERT_OK(delete F);
}

void TestIndBox::testCall3() {
    const size_t n = 30;
    Matrix lb(n, 1);
    Matrix ub(n, 1);
    for (int i = 0; i < n - 1; i++) {
        lb[i] = -i - 1;
        ub[i] = i + 1;
    }
    lb[n - 1] = -INFINITY;
    ub[n - 1] = INFINITY;
    Matrix x(n, 1);

    x[n - 1] = 1e9;

    Function *F = new IndBox(lb, ub);

    double f;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, f));

    _ASSERT_EQ(1.0, f);

    lb[n - 1] = 0.0;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, f));
    _ASSERT_EQ(1.0, f);

    x[n - 1] = -0.1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, f));
    _ASSERT(isinf(f));

    _ASSERT_OK(delete F);
}


