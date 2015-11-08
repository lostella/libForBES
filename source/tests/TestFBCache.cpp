#include "FBCache.h"
#include "FBProblem.h"

#include "TestFBCache.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TestFBCache);

TestFBCache::TestFBCache() {
}

TestFBCache::~TestFBCache() {
}

void TestFBCache::setUp() {
}

void TestFBCache::tearDown() {
}

void TestFBCache::test_boxqp_small() {
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
	// reference results 2
	double ref_fx2 = -6.000000000000000;
	double ref_y2_g01[] = {-0.50, -0.80, -1.30, -1.60};
	double ref_z2_g01[] = {-0.50, -0.80, -1.00, -1.00};
	double ref_y2_g02[] = {+0.00, -0.60, -1.60, -2.20};
	double ref_z2_g02[] = {+0.00, -0.60, -1.00, -1.00};
	double ref_FBEx2_g01 =  -7.450000000000000;
	double ref_FBEx2_g02 =  -8.900000000000000;

	Matrix Q = Matrix(n, n, data_Q);
	Matrix q = Matrix(n, 1, data_q);
	Quadratic f = Quadratic(Q, q);
	IndBox g = IndBox(lb, ub);
	FBProblem prob = FBProblem(f, g);
	
	FBCache * cache;
	Matrix * x, * y, * z;
	double fx, FBEx;
	
	// test FB operations starting from x1
	x = new Matrix(n, 1, data_x1);
	cache = new FBCache(prob, *x, 1.0);
	fx = cache->get_eval_f();
	y = cache->get_forward_step(gamma1);
	z = cache->get_forward_backward_step(gamma1);
	FBEx = cache->get_eval_FBE(gamma1);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx1, fx, 1e-14);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y1_g01[i], y->get(i, 0), 1e-14);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z1_g01[i], z->get(i, 0), 1e-14);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx1_g01, FBEx, 1e-14);
	
	// change gamma
	y = cache->get_forward_step(gamma2);
	z = cache->get_forward_backward_step(gamma2);
	FBEx = cache->get_eval_FBE(gamma2);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx1, fx, 1e-14);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y1_g02[i], y->get(i, 0), 1e-14);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z1_g02[i], z->get(i, 0), 1e-14);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx1_g02, FBEx, 1e-14);
	
	delete x;
	
	// change point x
	x = new Matrix(n, 1, data_x2);
	cache->set_x(*x);
	fx = cache->get_eval_f();
	y = cache->get_forward_step(gamma1);
	z = cache->get_forward_backward_step(gamma1);
	FBEx = cache->get_eval_FBE(gamma1);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx2, fx, 1e-14);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y2_g01[i], y->get(i, 0), 1e-14);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z2_g01[i], z->get(i, 0), 1e-14);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx2_g01, FBEx, 1e-14);
	
	// change gamma
	y = cache->get_forward_step(gamma2);
	z = cache->get_forward_backward_step(gamma2);
	FBEx = cache->get_eval_FBE(gamma2);
	
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx2, fx, 1e-14);
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_y2_g02[i], y->get(i, 0), 1e-14);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z2_g02[i], z->get(i, 0), 1e-14);
	}
	CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx2_g02, FBEx, 1e-14);
	
	delete x;
	delete cache;
}

void TestFBCache::test_l1logreg_small() {

}
