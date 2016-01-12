/*
 * File:   TestIndBox.cpp
 * Author: Pantelis Sopasakis
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

    Matrix x(2, 1);
    x[0] = -1.0;
    x[1] = 4.0;

    double fval = -1;
    _ASSERT(F->category().defines_f());
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_EQ(0.0, fval);

    x[1] = 5.0;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT(isinf(fval));

    x[0] = -0.5;
    x[1] = 3.1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_NOT(isinf(fval));
    _ASSERT_EQ(0.0, fval);

    x[0] = -1 - 1e-9;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_NEQ(0.0, fval);
    _ASSERT(isinf(fval));


    _ASSERT_OK(delete F);

}

void TestIndBox::testCall2() {
    const long n = 10;
    Matrix lb(n, 1);
    Matrix ub(n, 1);
    for (long i = 0; i < n; i++) {
        lb[i] = -i - 1;
        ub[i] = i + 1;
    }

    Function *F;
    F = new IndBox(lb, ub);

    Matrix x(n, 1);

    double fval = -1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_EQ(0.0, fval);

    x = ub;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_EQ(0.0, fval);

    x = lb;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_EQ(0.0, fval);


    Matrix y(n, n);
    _ASSERT_EXCEPTION(F->call(y, fval), std::invalid_argument);

    lb = Matrix(n, n);    
    _ASSERT_EXCEPTION(F = new IndBox(lb, ub), std::invalid_argument);
    
    
    lb = Matrix(n, 1);
    ub = Matrix(n, n);
   
    _ASSERT_EXCEPTION(F = new IndBox(lb, ub), std::invalid_argument);
    delete F;
}

void TestIndBox::testCall3() {
    const long n = 30;
    Matrix lb(n, 1);
    Matrix ub(n, 1);
    for (long i = 0; i < n - 1; i++) {
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

    _ASSERT_EQ(0.0, f);

    lb[n - 1] = 0.0;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, f));
    _ASSERT_EQ(0.0, f);

    x[n - 1] = -0.1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, f));
    _ASSERT(isinf(f));

    _ASSERT_OK(delete F);
}

void TestIndBox::testCallProx() {
    const double tol = 1e-10;
    double lb = -1.0;
    double ub = 4.0;
    Function *F = new IndBox(lb, ub);

    Matrix x(2, 1);
    x[0] = lb - 1.0;
    x[1] = (lb + ub) / 2.0;

    double gamma = 1.5;

    Matrix prox(2, 1);
    double fprox;

    _ASSERT(F->category().defines_prox());
    int status = F->callProx(x, gamma, prox, fprox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NUM_EQ(0.0, fprox, tol);

    _ASSERT_NUM_EQ(lb, prox.get(0, 0), tol);
    _ASSERT_NUM_EQ(x.get(1, 0), prox.get(1, 0), tol);

    status = F->callProx(x, gamma, prox);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    _ASSERT_NUM_EQ(lb, prox.get(0, 0), tol);
    _ASSERT_NUM_EQ(x.get(1, 0), prox.get(1, 0), tol);

    _ASSERT_OK(delete F);
}

void TestIndBox::testCategory() {
    double lb = 1.0;
    double ub = 2.0;
    Function * ind_box = new IndBox(lb, ub);
    FunctionOntologicalClass cat = ind_box -> category();
    _ASSERT(cat.defines_conjugate());
    _ASSERT(cat.defines_prox());
    _ASSERT(cat.defines_f());
    _ASSERT_NOT(cat.defines_conjugate_grad());
    _ASSERT_NOT(cat.defines_grad());
    _ASSERT_NOT(cat.defines_hessian());
    _ASSERT_NOT(cat.defines_hessian_conj());
    delete ind_box;
}

void TestIndBox::testCallConj() {
    double lb = -1.5;
    double ub = 2.0;
    Function * ind_box = new IndBox(lb, ub);
    
    Matrix x(5,1);
    x[1] = 10;
    x[2] = -0.5;
    x[3] = -1.0;
    x[4] = -1.5;
    double f_star = 0.0;
    ind_box->callConj(x, f_star);    
    _ASSERT_NUM_EQ(24.5, f_star, 1e-8);
    delete ind_box;
}


