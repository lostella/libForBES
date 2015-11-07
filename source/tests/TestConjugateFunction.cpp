/*
 * File:   TestConjugateFunction.cpp
 * Author: chung
 *
 * Created on Nov 7, 2015, 4:10:56 PM
 */

#include "TestConjugateFunction.h"
#include "ConjugateFunction.h"


#define FORBES_TEST_UTILS
#include "ForBES.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TestConjugateFunction);

TestConjugateFunction::TestConjugateFunction() {
}

TestConjugateFunction::~TestConjugateFunction() {
}

void TestConjugateFunction::setUp() {
}

void TestConjugateFunction::tearDown() {
}

void TestConjugateFunction::testCall() {
    Function *f = new QuadraticLoss(1.55628554);
    Function *f_conj = new ConjugateFunction(*f);

    Matrix x = MatrixFactory::MakeRandomMatrix(10, 1, 3.0, 5.0);

    double f_star;
    double f_conj_val;

    int status;

    _ASSERT(f->category().defines_f());

    status = f->callConj(x, f_star);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    status = f_conj->call(x, f_conj_val);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    _ASSERT_NUM_EQ(f_star, f_conj_val, 1e-8);

    delete f;
    delete f_conj;
}

void TestConjugateFunction::testCall2() {

    int n = 8;
    int s = 4;
    Matrix Q = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix A = MatrixFactory::MakeRandomMatrix(s, n, 0.0, -5.0, Matrix::MATRIX_DENSE);

    Matrix q = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix b = MatrixFactory::MakeRandomMatrix(s, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);

    QuadOverAffine * qoa;
    qoa = new QuadOverAffine(Q, q, A, b);
    _ASSERT_NEQ(NULL, qoa);

    Matrix y = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    double fstar = 0.0;
    Matrix grad;
    int status = qoa->callConj(y, fstar, grad);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_NOT(std::abs(fstar) < 1e-7);

    Function * f_conjugate = new ConjugateFunction(*qoa);
    double fstar2 = 0.0;
    Matrix grad2;
    f_conjugate -> call(y, fstar2, grad2);

    _ASSERT_NUM_EQ(fstar, fstar2, 1e-8);
    _ASSERT_EQ(grad, grad2);

    _ASSERT_OK(delete qoa);
    _ASSERT_OK(delete f_conjugate);


}

void TestConjugateFunction::testCallConj() {
    size_t n = 10;
    size_t nnz_Q = 20;
    Matrix Qsp = MatrixFactory::MakeRandomSparse(n, n, nnz_Q, 0.0, 1.0);

    Function *F = new Quadratic(Qsp);

    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 3.0, 1.5, Matrix::MATRIX_DENSE);
    double fval = -1.0;
    double fval2 = -1.0;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT(fval > 0);

    double f_exp = Qsp.quad(x);
    const double tol = 1e-10;
    _ASSERT_NUM_EQ(f_exp, fval, tol);


    Function * F_conj = new ConjugateFunction(*F);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F_conj->callConj(x, fval2));
    _ASSERT(fval2 > 0);
    _ASSERT_NUM_EQ(fval, fval2, tol);

    _ASSERT_OK(delete F);
    _ASSERT_OK(delete F_conj);
}

void TestConjugateFunction::testCallConj2() {

}

void TestConjugateFunction::testCallProx() {

}

void TestConjugateFunction::testCategory() {

}


