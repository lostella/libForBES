#include "FBCache.h"
#include "FBProblem.h"

#include "TestFBCache.h"

#define DOUBLES_EQUAL_DELTA 1e-8

CPPUNIT_TEST_SUITE_REGISTRATION(TestFBCache);

TestFBCache::TestFBCache() {
}

TestFBCache::~TestFBCache() {
}

void TestFBCache::setUp() {
}

void TestFBCache::tearDown() {
}

void TestFBCache::testBoxQP_small() {
	size_t n = 4;
	// problem data
	double data_Q[] = {
		7, 2, -2, -1,
		2, 3, 0, -1,
		-2, 0, 3, -1,
		-1, -1, -1, 1
	};
	double data_q[] = {
		1, 2, 3, 4
	};
	double gamma1 = 0.1;
	double gamma2 = 0.2;
	double lb = -1;
	double ub = +1;
	// starting points
	double data_x1[] = {+0.5, +1.2, -0.7, -1.1};
	double data_x2[] = {-1.0, -1.0, -1.0, -1.0};
	// reference results 1
	double ref_fx1 = +3.774999999999999;
	double ref_y1_g01[] = {-0.44, +0.43, -0.80, -1.29};
	double ref_z1_g01[] = {-0.44, +0.43, -0.80, -1.00};
	double ref_y1_g02[] = {-1.38, -0.34, -0.90, -1.48};
	double ref_z1_g02[] = {-1.00, -0.34, -0.90, -1.00};
	double ref_FBEx1_g01 =  -3.417500000000000;
	double ref_FBEx1_g02 = -10.514000000000003;
	/*
	 * WARNING: The following two arrays are never used
         */
	double ref_gradFBEx1_g01[] = {1.379999999999999, 3.409999999999999, 2.480000000000000, 0.909999999999999};
	double ref_gradFBEx1_g02[] = {-5.779999999999999, -0.020000000000000, 3.300000000000000, 2.840000000000000};

	// reference results 2
	double ref_fx2 = -6.000000000000000;
	double ref_y2_g01[] = {-0.50, -0.80, -1.30, -1.60};
	double ref_z2_g01[] = {-0.50, -0.80, -1.00, -1.00};
	double ref_y2_g02[] = {+0.00, -0.60, -1.60, -2.20};
	double ref_z2_g02[] = {+0.00, -0.60, -1.00, -1.00};
	double ref_FBEx2_g01 =  -7.450000000000000;
	double ref_FBEx2_g02 =  -8.900000000000000;

	Matrix Q(n, n, data_Q);
	Matrix q(n, 1, data_q);
	Quadratic f(Q, q);
	IndBox g(lb, ub);
	FBProblem prob(f, g);
	
	FBCache * cache;
	Matrix * x, * y, * z, * gradFBEx;
	double fx, FBEx;
	
	// test FB operations starting from x1
	x = new Matrix(n, 1, data_x1);
	cache = new FBCache(prob, *x, 1.0);
	fx = cache->get_eval_f();
	y = cache->get_forward_step(gamma1);
	z = cache->get_forward_backward_step(gamma1);
	FBEx = cache->get_eval_FBE(gamma1);
	gradFBEx = cache->get_grad_FBE(gamma1);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx1, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y1_g01[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z1_g01[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
		// CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_gradFBEx1_g01[i], gradFBEx->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx1_g01, FBEx, DOUBLES_EQUAL_DELTA);
	
	// change gamma
	y = cache->get_forward_step(gamma2);
	z = cache->get_forward_backward_step(gamma2);
	FBEx = cache->get_eval_FBE(gamma2);
	gradFBEx = cache->get_grad_FBE(gamma2);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx1, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y1_g02[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z1_g02[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
		// CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_gradFBEx1_g02[i], gradFBEx->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx1_g02, FBEx, DOUBLES_EQUAL_DELTA);
	
	delete x;
	
	// change point x
	x = new Matrix(n, 1, data_x2);
	cache->set_point(*x);
	fx = cache->get_eval_f();
	y = cache->get_forward_step(gamma1);
	z = cache->get_forward_backward_step(gamma1);
	FBEx = cache->get_eval_FBE(gamma1);
	gradFBEx = cache->get_grad_FBE(gamma1);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx2, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y2_g01[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z2_g01[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx2_g01, FBEx, DOUBLES_EQUAL_DELTA);
	
	// change gamma
	y = cache->get_forward_step(gamma2);
	z = cache->get_forward_backward_step(gamma2);
	FBEx = cache->get_eval_FBE(gamma2);
	gradFBEx = cache->get_grad_FBE(gamma2);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx2, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y2_g02[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z2_g02[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx2_g02, FBEx, DOUBLES_EQUAL_DELTA);
	
	delete x;
	delete cache;
}

void TestFBCache::testSparseLeastSquares_small() {
	size_t n = 10;
	size_t m = 5;
	// problem data
	double data_A[] = {
  		 4,    -5,     0,    -3,     1,
		-4,     2,     3,     8,    -1,
	   -11,    -5,     6,    -6,     4,
		 0,     7,   -10,    -1,    -7,
		14,     4,     6,    -6,    -3,
		-2,     5,    -2,     3,   -11,
		-2,    -5,    -8,     2,     1,
		 0,    -7,     5,     1,    -2,
		 0,    -2,    -9,    -2,    -5,
		-5,    -6,    -3,   -11,     4
	};
	double data_minusb[] = {1, 4, -6, 2, 3};
	double gamma1 = 0.0017;
	double gamma2 = 0.003;
	// starting points
	double data_x1[] = {-0.14, -0.24, -0.15, 0.03, 0.03, 0.04, -0.02, 0.01, -0.05, 0.08};
	double data_x2[] = {-0.14, 0.04,  -0.09, -0.04,0.05, 0.10, -0.12, 0.12, 0.06, -0.01};
	// reference results 1
	double ref_fx1 = 47.978750000000005;
	double ref_y1_g1[] = {-0.118172, -0.202073, 0.015444, -0.128984, 0.015125, 0.021640, -0.072156, 0.141121, -0.124545, 0.114323};
	double ref_z1_g1[] = {-0.116472, -0.200373, 0.013744, -0.127284, 0.013425, 0.019940, -0.070456, 0.139421, -0.122845, 0.112623};
	double ref_y1_g2[] = {-0.10148,-0.17307,0.14196,-0.25056,0.00375,0.00760,-0.11204, 0.24139,-0.18155,0.14057};
	double ref_z1_g2[] = {-0.09848,-0.17007,0.13896,-0.24756,0.00075,0.00460,-0.10904, 0.23839,-0.17855,0.13757};
	double ref_FBEx1_g1 =  24.874158335000004;
	double ref_FBEx1_g2 =  6.877117650000010;
	// reference results 2
	double ref_fx2 = 33.691600000000001;
	double ref_y2_g1[] = {-0.093607, 0.018614, 0.070837, -0.171631, 0.050833, 0.043475, -0.145942, 0.222323, 0.021240, 0.091524};
	double ref_z2_g1[] = {  -0.091907, 0.016914, 0.069137, -0.169931, 0.049133, 0.041775, -0.144242, 0.220623, 0.019540, 0.089824};
	double ref_y2_g2[] = {-0.05813,0.00226,0.19383,-0.27229,0.05147,0.00025,-0.16578,0.30057,-0.00840,0.16916};
	double ref_z2_g2[] = {-0.05513,0,0.19083,-0.26929,0.04847,0,-0.16278,0.29757,-0.00540,0.16616};
	double ref_FBEx2_g1 =  13.450436129999996;
	double ref_FBEx2_g2 =  -2.445831616666666;

	Matrix A(m, n, data_A);
	Matrix minusb(m, 1, data_minusb);
	MatrixOperator OpA(A);
	QuadraticLoss f();
	Norm1 g();
	FBProblem prob(f, OpA, minusb, g);
	
	FBCache * cache;
	Matrix * x, * y, * z, * gradFBEx;
	double fx, FBEx;
	
	// test FB operations starting from x1
	x = new Matrix(n, 1, data_x1);
	cache = new FBCache(prob, *x, 1.0);
	fx = cache->get_eval_f();
	y = cache->get_forward_step(gamma1);
	z = cache->get_forward_backward_step(gamma1);
	FBEx = cache->get_eval_FBE(gamma1);
	gradFBEx = cache->get_grad_FBE(gamma1);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx1, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y1_g1[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z1_g1[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx1_g1, FBEx, DOUBLES_EQUAL_DELTA);
	
	// change gamma
	y = cache->get_forward_step(gamma2);
	z = cache->get_forward_backward_step(gamma2);
	FBEx = cache->get_eval_FBE(gamma2);
	gradFBEx = cache->get_grad_FBE(gamma2);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx1, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y1_g2[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z1_g2[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx1_g2, FBEx, DOUBLES_EQUAL_DELTA);
	
	delete x;
	
	// test FB operations starting from x2
	x = new Matrix(n, 1, data_x2);
	cache->set_point(*x);
	fx = cache->get_eval_f();
	y = cache->get_forward_step(gamma1);
	z = cache->get_forward_backward_step(gamma1);
	FBEx = cache->get_eval_FBE(gamma1);
	gradFBEx = cache->get_grad_FBE(gamma1);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx2, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y2_g1[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z2_g1[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx2_g1, FBEx, DOUBLES_EQUAL_DELTA);
	
	// change gamma
	y = cache->get_forward_step(gamma2);
	z = cache->get_forward_backward_step(gamma2);
	FBEx = cache->get_eval_FBE(gamma2);
	gradFBEx = cache->get_grad_FBE(gamma2);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx2, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y2_g2[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z2_g2[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx2_g2, FBEx, DOUBLES_EQUAL_DELTA);
	
	delete x;	
	delete cache;
}

void TestFBCache::testSparseLogReg_small() {
	size_t n = 5;
	size_t m = 4;
	// problem data
	double data_A[] = {
		1, 2, -1, -1,
		-2, -1, 0, -1,
		3, 0, 4, -1,
		-4, -1, -3, 1,
		5, 3, 2, 3
	};
	double data_minusb[] = {
		-1, 1, -1, 1
	};
	double gamma1 = 0.1;
	double gamma2 = 0.05;
	// starting points
	double data_x1[] = {0, 0, 0, 0, 0};
	// reference results
	double ref_fx1 = 3.253046750072892;
	double ref_y1_g1[] = {0.026894142137000, -0.200000000000000, 0.484846862904004, -0.511741005041004, 0.673105857863000};
	double ref_z1_g1[] = {0, -0.100000000000000, 0.384846862904004, -0.411741005041004, 0.573105857863001};
	double ref_FBEx1_g1 = -0.027393687107681;
	double ref_y1_g2[] = {0.013447071068500, -0.100000000000000, 0.242423431452002, -0.255870502520502, 0.336552928931500};
	double ref_z1_g2[] = {0, -0.050000000000000, 0.192423431452002, -0.205870502520502, 0.286552928931500};
	double ref_FBEx1_g2 = 1.612826531482605;
	
	Matrix A(m, n, data_A);
	Matrix minusb(m, 1, data_minusb);
	Matrix xstar;
	LogLogisticLoss f(1.0);
	MatrixOperator OpA(A);
	Norm1 g(1.0);
	FBProblem prob(f, OpA, minusb, g);

	FBCache * cache;
	Matrix * x, * y, * z, * gradFBEx;
	double fx, FBEx;
	
	// test FB operations starting from x1
	x = new Matrix(n, 1, data_x1);
	cache = new FBCache(prob, *x, 1.0);
	fx = cache->get_eval_f();
	y = cache->get_forward_step(gamma1);
	z = cache->get_forward_backward_step(gamma1);
	FBEx = cache->get_eval_FBE(gamma1);
	gradFBEx = cache->get_grad_FBE(gamma1);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx1, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y1_g1[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z1_g1[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx1_g1, FBEx, DOUBLES_EQUAL_DELTA);
	
	// change gamma
	y = cache->get_forward_step(gamma2);
	z = cache->get_forward_backward_step(gamma2);
	FBEx = cache->get_eval_FBE(gamma2);
	gradFBEx = cache->get_grad_FBE(gamma2);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx1, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y1_g2[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z1_g2[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx1_g2, FBEx, DOUBLES_EQUAL_DELTA);
	
	delete x;

	delete cache;
}

void TestFBCache::testLogLossPlusL1_small() {
	size_t n = 5;
	size_t m = 5;
	double gamma1 = 0.1;
	double gamma2 = 0.05;
	// starting points
	double data_x1[] = {1, -2, 3, -4, 5};
	// reference results
	double ref_fx1 = 6.513642326541865;
	double ref_y1_g1[] = {1.026894142136999, -1.911920292202212, 3.004742587317756, -3.901798620996209, 5.000669285092428};
	double ref_z1_g1[] = {0.926894142136999, -1.811920292202212, 2.904742587317756, -3.801798620996209, 4.900669285092428};
	double ref_FBEx1_g1 = 21.018928350676088;
	double ref_gradFBEx1_g1[] = {0.716685094584284, -1.861049915114549, 0.948270715102844, -1.978513017309498, 0.992646792853859};
	double ref_y1_g2[] = {1.013447071068500, -1.955960146101106, 3.002371293658878, -3.950899310498105, 5.000334642546214};
	double ref_z1_g2[] = {0.963447071068500, -1.905960146101106, 2.952371293658878, -3.900899310498105, 4.950334642546214};
	double ref_FBEx1_g2 = 21.266285338608977;
	double ref_gradFBEx1_g2[] = {0.723871836607145, -1.870923496546218, 0.950422420962636, -1.980263403673695, 0.992976970964787};

	Function * f = new LogLogisticLoss(1.0);
	Function * g = new Norm1(1.0);
	FBProblem * prob = new FBProblem(*f, *g);

	FBCache * cache;
	Matrix * x, * y, * z, * gradFBEx;
	double fx, FBEx;
	
	// test FB operations starting from x1
	x = new Matrix(n, 1, data_x1);
	cache = new FBCache(*prob, *x, 1.0);
	fx = cache->get_eval_f();
	y = cache->get_forward_step(gamma1);
	z = cache->get_forward_backward_step(gamma1);
	FBEx = cache->get_eval_FBE(gamma1);
	gradFBEx = cache->get_grad_FBE(gamma1);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx1, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y1_g1[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z1_g1[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_gradFBEx1_g1[i], gradFBEx->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx1_g1, FBEx, DOUBLES_EQUAL_DELTA);
	
	// change gamma
	y = cache->get_forward_step(gamma2);
	z = cache->get_forward_backward_step(gamma2);
	FBEx = cache->get_eval_FBE(gamma2);
	gradFBEx = cache->get_grad_FBE(gamma2);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx1, fx, DOUBLES_EQUAL_DELTA);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y1_g2[i], y->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z1_g2[i], z->get(i, 0), DOUBLES_EQUAL_DELTA);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_gradFBEx1_g2[i], gradFBEx->get(i, 0), DOUBLES_EQUAL_DELTA);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx1_g2, FBEx, DOUBLES_EQUAL_DELTA);
	
	delete x;
	delete cache;

	delete prob;
	delete f;
	delete g;
}

