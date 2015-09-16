/*
 * File:   TestQuadratic.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Jul 9, 2015, 4:14:39 AM
 */

#include "TestQuadratic.h"


const static double MAT1[16] = {
    7, 2, -2, -1,
    2, 3, 0, -1,
    -2, 0, 3, -1,
    -1, -1, -1, 1
};
const static double MAT2[16] = {
    16.0, 2.0, 3.0, 13.0,
    5.0, 11.0, 10.0, 8.0,
    9.0, 7.0, 6.0, 12.0,
    4.0, 14.0, 15.0, 1.0
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestQuadratic);

TestQuadratic::TestQuadratic() {

}

TestQuadratic::~TestQuadratic() {
}

void TestQuadratic::setUp() {
}

void TestQuadratic::tearDown() {
}

void TestQuadratic::testQuadratic() {
    const double * Qdata;
    Qdata = MAT2;
    double xdata[4] = {-1.0, 1.0, 1.0, 1.0};

    Matrix Q = Matrix(4, 4, Qdata);
    Matrix x = Matrix(4, 1, xdata);
    Function *quad = NULL;
    _ASSERT_OK(quad = new Quadratic(Q));
    _ASSERT_NEQ(quad, NULL);
    
    double f;
    int info;
    _ASSERT_OK(info = quad -> call(x, f));
    _ASSERT_EQ(ForBESUtils::STATUS_OK, info);

    if (quad != NULL) {
        delete quad;
    }
}

void TestQuadratic::testQuadratic2() {
    /* Test the empty constructor */
    size_t n = 8;
    const double tol = 1e-7;
    Function *F = new Quadratic(); /* Q = I, q = 0. */
    _ASSERT_EQ(Function::CAT_QUADRATIC, F->category());

    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    x.set(0, 0, 666);
    Matrix grad;

    double fval = 0.0;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval, grad));

    const double fval_expected = (x * x).get(0, 0) / 2;
    _ASSERT_NUM_EQ(fval_expected, fval, tol);

    Matrix q(n, 1);
    for (size_t i = 0; i < n; ++i) {
        _ASSERT_OK(q.set(i, 0, i + 1));
    }
    _ASSERT_OK(static_cast<Quadratic*> (F)->setq(q));

    double fval_prev = fval;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval, grad));
    _ASSERT_NUM_EQ(fval_prev + (q * x).get(0, 0), fval, 1e-6);

    _ASSERT_OK(delete F);
}

void TestQuadratic::testQuadratic3() {
    size_t n = 10;
    Matrix Q(n, n, Matrix::MATRIX_DENSE);
    for (size_t i = 0; i < n - 1; ++i) {
        Q.set(i, i, 1.5);
    }
    Quadratic *F = new Quadratic(Q);
    Matrix x(n, 1);
    double fval;
    int status = F->callConj(x, fval);
    _ASSERT_EQ(ForBESUtils::STATUS_NUMERICAL_PROBLEMS, status);

    Q.set(n - 1, n - 1, 1.5);
    Q.set(1, 0, 0.1);
    Q.set(0, 1, 0.1);

    F->setQ(Q);
    _ASSERT_OK(status = F->callConj(x, fval));
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    //
    //    Matrix q(n, 1);
    //    for (size_t i = 0; i < n; ++i) {
    //        q.set(i, 0, i + 1);
    //        x.set(i, 0, 2 * i + 1);
    //    }
    //
    //
    //    Matrix grad;
    //    static_cast<Quadratic*> (F)->setq(q);
    //
    //    _ASSERT_OK(status = F->callConj(x, fval, grad));
    //    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    //
    //    const double fval_exp = 190.002976190476;
    //    const double tol = 1e-8;
    //
    //    _ASSERT_NUM_EQ(fval_exp, fval, tol);
    //
    //    double grad0_exp = -0.0446428571428572;
    //    _ASSERT_NUM_EQ(grad0_exp, grad.get(0, 0), tol);
    //
    //    std::cout << fval;
    //    _ASSERT_OK(delete F);
}

void TestQuadratic::testCallProx() {
    Function *F = new Quadratic();
    Matrix x;
    Matrix prox;
    _ASSERT_EQ(ForBESUtils::STATUS_UNDEFINED_FUNCTION, F->callProx(x, 0.5, prox));
    _ASSERT_OK(delete F);
}

void TestQuadratic::testCall() {
    const double * Qdata;
    Qdata = MAT2;
    double qdata[4] = {2.0, 3.0, 4.0, 5.0};
    double xdata[4] = {-1.0, 1.0, 1.0, 1.0};

    Matrix Q = Matrix(4, 4, Qdata);
    Matrix q = Matrix(4, 1, qdata);
    Matrix x = Matrix(4, 1, xdata);

    Quadratic quadratic(Q, q);
    double f = -999.0;
    int status = quadratic.call(x, f);

    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_EQ(42.0, f);

    /* Second part */
    Quadratic quadratic2(Q);
    status = quadratic2.call(x, f);
    _ASSERT_EQ(0, status);
    _ASSERT_EQ(32.0, f);
}

void TestQuadratic::testCallWithGradient() {
    const size_t n = 4;
    const double * Qdata;
    Qdata = MAT1;
    double xdata[4] = {-1.0, 1.0, 1.0, 1.0};
    const double expected[4] = {-8, 0, 4, 0};

    Matrix Q = Matrix(n, n, Qdata);
    Matrix x = Matrix(n, 1, xdata);

    Function * quad = new Quadratic(Q);
    double f = -999.0f;
    Matrix grad;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, quad -> call(x, f, grad));

    for (size_t i = 0; i < n; i++) {
        _ASSERT_EQ(expected[i], grad[i]);
    }

    delete quad;

}

void TestQuadratic::testCallConj() {
    const size_t n = 4;
    const double * Qdata;
    Qdata = MAT1;

    double qdata[4] = {2.0, 3.0, 4.0, 5.0};
    double xdata[4] = {-1.0, 1.0, 1.0, 1.0};

    Matrix Q = Matrix(n, n, Qdata);
    Matrix q = Matrix(n, 1, qdata);
    Matrix x = Matrix(n, 1, xdata);

    Quadratic quadratic(Q, q);
    double fstar;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, quadratic.callConj(x, fstar));


    double expected = 421.0;
    const double rel_tol = 1e-5;
    _ASSERT(std::fabs(expected - fstar) / expected < rel_tol);

    for (size_t i = 0; i < n; i++) {
        x[i] = 10 + 2 * i;
    }

    _ASSERT_EQ(ForBESUtils::STATUS_OK, quadratic.callConj(x, fstar));
    expected = 3722;
    _ASSERT(std::fabs(expected - fstar) / expected < rel_tol);

    Matrix grad;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, quadratic.callConj(x, fstar, grad));

    Matrix grad_expected(n, 1);
    grad_expected.set(0, 0, 45.5);
    grad_expected.set(1, 0, 35.5);
    grad_expected.set(2, 0, 96.5);
    grad_expected.set(3, 0, 188.5);

    const double tol = 1e-5;
    for (size_t i = 0; i < n; i++) {
        _ASSERT_NUM_EQ(grad_expected.get(i, 0), grad.get(i, 0), tol);
    }

}

void TestQuadratic::testCallConj2() {
    const size_t n = 4;

    double qdata[4] = {2.0, 3.0, 4.0, 5.0};
    double xdata[4] = {-1.0, 1.0, 1.0, 1.0};

    Quadratic *F = new Quadratic();
    Matrix q(n, 1, qdata);
    Matrix x(n, 1, xdata);

    _ASSERT_OK(F->setq(q));
    double fstar;
    double fstar_exp = 38.0;
    const double tol = 1e-6;

    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->callConj(x, fstar));
    _ASSERT_NUM_EQ(fstar_exp, fstar, tol);
    _ASSERT_OK(delete F);

    double alpha = 2.3;
    Matrix Eye = MatrixFactory::MakeIdentity(n, alpha); /* Here Eye is a diagonal matrix */
    Function *F2 = new Quadratic(Eye, q);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F2->callConj(x, fstar));
    _ASSERT_NUM_EQ(fstar_exp / alpha, fstar, tol);
    _ASSERT_OK(delete F2);
}

void TestQuadratic::testCategory() {
    Quadratic quadratic;
    int cat = quadratic.category();
    _ASSERT_EQ(Function::CAT_QUADRATIC, cat);
}

void TestQuadratic::testCallDiagonalMatrix() {
    const size_t n = 4;
    Matrix Q(n, n, Matrix::MATRIX_DIAGONAL);
    for (size_t i = 0; i < n; i++) {
        Q[i] = i + 2.0f;
    }

    Function *f = new Quadratic(Q);

    Matrix x(n, 1);
    for (size_t i = 0; i < n; i++) {
        x[i] = 2 * i + 1.0f;
    }

    double val;
    f -> call(x, val);

    double tol = 1e-8;
    _ASSERT_NUM_EQ(187.0, val, tol);
    _ASSERT_OK(delete f);
}

void TestQuadratic::testCallSparse() {
    /*
     * Q is sparse, x is dense
     */
    size_t n = 10;
    size_t nnz_Q = 20;
    Matrix Qsp = MatrixFactory::MakeRandomSparse(n, n, nnz_Q, 0.0, 1.0);

    Function *F = new Quadratic(Qsp);

    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 3.0, 1.5, Matrix::MATRIX_DENSE);
    double fval = -1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT(fval > 0);

    double f_exp = Qsp.quad(x);
    const double tol = 1e-8;
    _ASSERT_NUM_EQ(f_exp, fval, tol);

    _ASSERT_OK(delete F);

}

void TestQuadratic::testCallSparse2() {
    /*
     * BOTH Q and x are sparse
     */
    size_t n = 10;
    size_t nnz_Q = 20;
    size_t nnz_x = 6;
    Matrix Qsp = MatrixFactory::MakeRandomSparse(n, n, nnz_Q, 0.0, 1.0);

    Function *F = new Quadratic(Qsp);

    Matrix x = MatrixFactory::MakeRandomSparse(n, 1, nnz_x, 0.0, 1.0);
    double fval = -1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT(fval > 0);

    double f_exp = Qsp.quad(x);
    const double tol = 1e-8;
    _ASSERT_NUM_EQ(f_exp, fval, tol);

    _ASSERT_OK(delete F);
}

void TestQuadratic::testCallConjSparse() {
    size_t n = 5;
    size_t nnz_Q = 2 * n - 1;
    double fstar = -1;
    const double tol = 1e-8;
    const double fstar_exp = 5.13142364727941;
    int status;
    Function *F;

    Matrix Qsp = MatrixFactory::MakeSparseSymmetric(n, nnz_Q);
    Matrix x(n, 1);

    for (size_t i = 0; i < n; i++) {
        Qsp.set(i, i, 10.0);
        x.set(i, 0, i + 1);
    }
    for (size_t i = 1; i < n; i++) { /* Set the LT part */
        Qsp.set(i, i - 1, 0.5);
    }

    F = new Quadratic(Qsp);
    _ASSERT_OK(status = F->callConj(x, fstar));
    _ASSERT_NEQ(-1, fstar);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

    _ASSERT_NUM_EQ(fstar_exp, fstar, tol);

    _ASSERT_OK(delete F);

}

void TestQuadratic::testCallSparse3() {
    const size_t n = 10;
    const double tol = 1e-7;
    Matrix Q = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 2.0, Matrix::MATRIX_SYMMETRIC);
    Matrix q = MatrixFactory::MakeSparse(n, 1, 0, Matrix::SPARSE_UNSYMMETRIC); /* q = [] */
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);

    Function *F = new Quadratic(Q, q);

    double f;
    int status = F->call(x, f);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    double f_exp = Q.quad(x);
    _ASSERT_NUM_EQ(f_exp, f, tol);

    q = MatrixFactory::MakeSparse(n, 1, 1, Matrix::SPARSE_UNSYMMETRIC);
    q.set(0, 0, 5.5);


    status = F->call(x, f);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    f_exp = Q.quad(x, q);
    _ASSERT_NUM_EQ(f_exp, f, tol);

    _ASSERT_OK(delete F);

}