#include "FBProblem.h"

#include "TestLasso.h"

#define DOUBLES_EQUAL_DELTA 1e-8

CPPUNIT_TEST_SUITE_REGISTRATION(TestLasso);

TestLasso::TestLasso() {
}

TestLasso::~TestLasso() {
}

void TestLasso::setUp() {
}

void TestLasso::tearDown() {
}

void TestLasso::runTest() {
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

	Matrix A(m, n, data_A);
	Matrix minusb(m, 1, data_minusb);
	MatrixOperator OpA(A);
	QuadraticLoss f;
	Norm1 g;

	size_t n_tests = 100;
	for (int i_loop=0; i_loop<n_tests; i_loop++) {
		FBProblem prob(f, OpA, minusb, g);
	}
}