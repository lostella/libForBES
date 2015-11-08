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

void TestFBCache::test_boxqp_fx_z() {
	size_t n1 = 4;
	double data_Q1[] = {
		7, 2, -2, -1,
		2, 3, 0, -1,
		-2, 0, 3, -1,
		-1, -1, -1, 1
	};
	double data_q1[] = {
		1, 2, 3, 4
	};
	double lb1 = -1;
    double ub1 = +1;
    double data_x1[] = {
		0.5, 1.2, -0.7, -1.1
	};
	double cost1 = 3.774999999999999;
	double data_z1[] = {
		-0.44, 0.43, -0.8, -1
	};

    Matrix x = Matrix(n1, 1, data_x1);
    
    Matrix Q = Matrix(n1, n1, data_Q1);
    Matrix q = Matrix(n1, 1, data_q1);

    Quadratic f = Quadratic(Q, q);
    IndBox g = IndBox(lb1, ub1);
    
    FBProblem prob = FBProblem(f, g);
    
    FBCache cache = FBCache(prob, x, 1.0);
    
    double fx = cache.get_eval_f();
    
    Matrix * z = cache.get_forward_backward_step(0.1);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(cost1, fx, 1e-14);
  
    for (int i=0; i < n1; i++) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(data_z1[i], z->get(i, 0), 1e-14);
    }
}
