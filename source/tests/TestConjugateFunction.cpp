/*
 * File:   TestConjugateFunction.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Nov 7, 2015, 4:10:56 PM
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

#include "TestConjugateFunction.h"
#include "ConjugateFunction.h"
#include <cmath>

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
    _ASSERT(f_conj->category().defines_f());
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
    _ASSERT(f_conjugate->category().defines_grad());
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
    _ASSERT(F_conj->category().defines_conjugate());
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F_conj->callConj(x, fval2));
    _ASSERT(fval2 > 0);
    _ASSERT_NUM_EQ(fval, fval2, tol);

    _ASSERT_OK(delete F);
    _ASSERT_OK(delete F_conj);
}

void TestConjugateFunction::testCallConj2() {
    const size_t n = 10;
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.5, 2.0);
    const double delta = 0.2;
    Function * huber = new HuberLoss(delta);

    double f;
    double f2;
    Matrix grad(n, 1);
    Matrix grad2(n, 1);

    _ASSERT(huber->category().defines_f());
    int status = huber->call(x, f, grad);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    Function * huber_conj = new ConjugateFunction(*huber);
    _ASSERT(huber_conj->category().defines_conjugate());
    _ASSERT(huber_conj->category().defines_conjugate_grad());
    status = huber_conj -> callConj(x, f2, grad2);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    const double tol = 1e-9;
    _ASSERT_NUM_EQ(f, f2, tol);
}

void TestConjugateFunction::testCallProx() {
    Function * elastic = new ElasticNet(2.5, 1.3);
    Function * elastic_conj = new ConjugateFunction(*elastic);

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

    Matrix prox_conj(n, 1);
    status = elastic_conj -> callProx(x, gamma, prox_conj);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    
    
    
    delete elastic;
    delete elastic_conj;
}

void TestConjugateFunction::testCategory() {

}


