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
	double gamma = 0.1;
	double lb = -1;
    double ub = +1;
    // starting points
    double data_x1[] = {+0.5, +1.2, -0.7, -1.1};
    double data_x2[] = {-1.0, -1.0, -1.0, -1.0};
	// reference results
	double ref_fx1 = +3.774999999999999;
	double ref_fx2 = -6.000000000000000;
	double ref_z1[] = {-0.44, +0.43, -0.80, -1.00};
	double ref_z2[] = {-0.50, -0.80, -1.00, -1.00};
	double ref_FBEx1 = -3.417500000000000;
	double ref_FBEx2 = -7.450000000000000;

    Matrix Q = Matrix(n, n, data_Q);
    Matrix q = Matrix(n, 1, data_q);
    Quadratic f = Quadratic(Q, q);
    IndBox g = IndBox(lb, ub);
    FBProblem prob = FBProblem(f, g);
    
    FBCache * cache;
    Matrix * x, * z;
    double fx, FBEx;
    
    // test FB operations starting from x1
    x = new Matrix(n, 1, data_x1);
    cache = new FBCache(prob, *x, 1.0);
    fx = cache->get_eval_f();
    z = cache->get_forward_backward_step(gamma);
    FBEx = cache->get_eval_FBE(gamma);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx1, fx, 1e-14);
    for (int i=0; i < n; i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z1[i], z->get(i, 0), 1e-14);
    }
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx1, FBEx, 1e-14);
    delete x;
    delete cache;
    
    // test FB operations starting from x2
    x = new Matrix(n, 1, data_x2);
    cache = new FBCache(prob, *x, 1.0);
    fx = cache->get_eval_f();
    z = cache->get_forward_backward_step(gamma);
    FBEx = cache->get_eval_FBE(gamma);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_fx2, fx, 1e-14);
    for (int i=0; i < n; i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_z2[i], z->get(i, 0), 1e-14);
    }
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_FBEx2, FBEx, 1e-14);
    delete x;
    delete cache;
}

void TestFBCache::test_l1logreg_small() {

}
