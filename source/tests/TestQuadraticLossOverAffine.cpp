/*
 * File:   TestQuadraticLossOverAffine.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Nov 5, 2015, 12:06:30 AM
 */

#include "TestQuadraticLossOverAffine.h"
#include "MatrixFactory.h"
#include "QuadraticLossOverAffine.h"
#include "ForBES.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestQuadraticLossOverAffine);

TestQuadraticLossOverAffine::TestQuadraticLossOverAffine() {
}

TestQuadraticLossOverAffine::~TestQuadraticLossOverAffine() {
}

void TestQuadraticLossOverAffine::setUp() {
}

void TestQuadraticLossOverAffine::tearDown() {
}

void TestQuadraticLossOverAffine::testMethod() {
    size_t n = 8;
    size_t m = 2;
    
    Matrix A = MatrixFactory::MakeRandomMatrix(m, n, 0.0, 1.0);
    Matrix b = MatrixFactory::MakeRandomMatrix(m, 1, 0.0, 1.0);
    
    Matrix w = MatrixFactory::MakeRandomMatrix(n, 1, 10.0, 1.0);
    Matrix p = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0);
    
    Function *fun = new QuadraticLossOverAffine(A, b, w, p);
    
    Matrix y = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0);
    
    double f_star;
    Matrix grad(n,1);
    
    int status = fun->callConj(y, f_star, grad);
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
        
    delete fun;
   
}



