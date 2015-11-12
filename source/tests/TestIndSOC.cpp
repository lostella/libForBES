/*
 * File:   TestIndSOC.cpp
 * Author: Lorenzo Stella
 *
 * Created on Sept 21, 2015
 */

#include <math.h>
#include <cmath>
#include "TestIndSOC.h"
#include "Matrix.h"
#include "MatrixFactory.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TestIndSOC);

TestIndSOC::TestIndSOC() {
}

TestIndSOC::~TestIndSOC() {
}

void TestIndSOC::setUp() {
}

void TestIndSOC::tearDown() {
}

void TestIndSOC::testCall() {
    Function * F = new IndSOC(5);

    Matrix x(5, 1);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    x[3] = 4.0;
    x[4] = 5.0;

    double fval = -1;
    _ASSERT(F->category().defines_f());
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT(isinf(fval));

    x[4] = 6.0;
    fval = -1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT_EQ(0.0, fval);

    x[4] = -6.0;
    fval = -1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, F->call(x, fval));
    _ASSERT(isinf(fval));

    _ASSERT_OK(delete F);
}

void TestIndSOC::testCallProx() {
    Function * F = new IndSOC(5);
    int status;
    int eqflag;
    double fval;
    

    Matrix x(5, 1);
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    x[3] = 4.0;
    x[4] = 5.0;
    
    Matrix ref(5, 1);
    ref[0] = 9.564354645876385e-01;
    ref[1] = 1.912870929175277e+00;
    ref[2] = 2.869306393762916e+00;
    ref[3] = 3.825741858350554e+00;
    ref[4] = 5.238612787525831e+00;
    
    Matrix y(5, 1);
    fval = -1;
    _ASSERT(F->category().defines_prox());
    status = F->callProx(x, 1.0, y, fval);
    eqflag = 1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_EQ(0.0, fval);
    for (int i=0; i<5; i++) {
    	if (abs(y[i]-ref[i]) >= 1e-14) {
    		eqflag = 0;
    		break;
    	}
	}
	_ASSERT_EQ(1, eqflag);

    x[4] = 6.0;
    fval = -1;
    status = F->callProx(x, 1.0, y, fval);
    eqflag = 1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_EQ(0.0, fval);
    for (int i=0; i<5; i++) {
    	if (x[i] != y[i]) {
    		eqflag = 0;
    		break;
    	}
	}
	_ASSERT_EQ(1, eqflag);

    x[4] = -6.0;
    fval = -1;
    status = F->callProx(x, 1.0, y, fval);
    eqflag = 1;
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_EQ(0.0, fval);
	for (int i=0; i<5; i++) {
    	if (y[i] != 0.0) {
    		eqflag = 0;
    		break;
    	}
	}
	_ASSERT_EQ(1, eqflag);
	
    _ASSERT_OK(delete F);
}


