#include "FBSplitting.h"
#include "FBProblem.h"

#include "TestFBSplitting.h"

#include <iostream>

CPPUNIT_TEST_SUITE_REGISTRATION(TestFBSplitting);

TestFBSplitting::TestFBSplitting() {
}

TestFBSplitting::~TestFBSplitting() {
}

void TestFBSplitting::setUp() {
}

void TestFBSplitting::tearDown() {
}

void TestFBSplitting::testBoxQP_small() {
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
	double ref_xstar[] = {-0.352941176470588, -0.764705882352941, -1.000000000000000, -1.000000000000000};

	Matrix Q = Matrix(n, n, data_Q);
	Matrix q = Matrix(n, 1, data_q);
	Matrix * x0;
	Matrix xstar;
	Quadratic f = Quadratic(Q, q);
	IndBox g = IndBox(lb, ub);
	FBProblem prob = FBProblem(f, g);
	FBSplitting * solver;
	
	// test FB operations starting from x1
	x0 = new Matrix(n, 1, data_x1);
	solver = new FBSplitting(prob, *x0, gamma);
	solver->setMaxIt(1000);
	solver->setTol(1e-14);
	solver->run();
	xstar = solver->getSolution();
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_xstar[i], xstar.get(i, 0), 1e-14);
	}
	delete x0;
	delete solver;
	
	// test FB operations starting from x2
	x0 = new Matrix(n, 1, data_x2);
	solver = new FBSplitting(prob, *x0, gamma);
	solver->setMaxIt(1000);
	solver->setTol(1e-14);
	solver->run();
	xstar = solver->getSolution();
	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_xstar[i], xstar.get(i, 0), 1e-14);
	}
	delete x0;
	delete solver;
}
