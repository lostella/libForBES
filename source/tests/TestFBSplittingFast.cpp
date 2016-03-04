#include "FBSplittingFast.h"
#include "FBProblem.h"
#include "FBStoppingRelative.h"
#include "MatrixFactory.h"
#include "MatrixOperator.h"
#include "TestFBSplittingFast.h"

// #include <iostream>

#define DOUBLES_EQUAL_DELTA 1e-4
#define MAXIT 1000
#define TOLERANCE 1e-6

CPPUNIT_TEST_SUITE_REGISTRATION(TestFBSplittingFast);

TestFBSplittingFast::TestFBSplittingFast() {
}

TestFBSplittingFast::~TestFBSplittingFast() {
}

void TestFBSplittingFast::setUp() {
}

void TestFBSplittingFast::tearDown() {
}

void TestFBSplittingFast::testBoxQP_small() {
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

	Matrix * Q = new Matrix(n, n, data_Q);
	Matrix * q = new Matrix(n, 1, data_q);
	Matrix * x0;
	Matrix xstar;
	Function * f = new Quadratic(*Q, *q);
	Function * g = new IndBox(lb, ub);
	FBProblem prob = FBProblem(*f, *g);
	FBStoppingRelative sc = FBStoppingRelative(TOLERANCE);
	FBSplitting * solver;
	
	// test FB operations starting from x1
	size_t repeat = 100;
	for (size_t r = 0; r < repeat; r++) {
		x0 = new Matrix(n, 1, data_x1);
		solver = new FBSplittingFast(prob, *x0, gamma, sc, MAXIT);
		solver->run();
		xstar = solver->getSolution();
		//cout << "*** iters (fast) : " << solver->getIt() << endl;
		_ASSERT(solver->getIt() < MAXIT);
		for (int i=0; i < n; i++) {
			CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_xstar[i], xstar.get(i, 0), DOUBLES_EQUAL_DELTA);
		}
		delete x0;
		delete solver;
	
		// test FB operations starting from x2
		x0 = new Matrix(n, 1, data_x2);
		solver = new FBSplittingFast(prob, *x0, gamma, sc, MAXIT);
		solver->run();
		xstar = solver->getSolution();
		//cout << "*** iters (fast) : " << solver->getIt() << endl;
		_ASSERT(solver->getIt() < MAXIT);
		for (int i=0; i < n; i++) {
			CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_xstar[i], xstar.get(i, 0), DOUBLES_EQUAL_DELTA);
		}
		delete x0;
		delete solver;
	}

	delete Q;
	delete q;
	delete f;
	delete g;
}

void TestFBSplittingFast::testLasso_small() {
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
		-1, -2, -3, -4
	};
	/*
	 * WARNING: data_w is not used anywhere...
	 */
	double data_w[] = {
		1, 1, 1, 1
	};
	double gamma = 0.01;
	// starting points
	double data_x1[] = {0, 0, 0, 0, 0};
	// reference results
	double ref_xstar[] = {-0.010238907849511, 0, 0, 0, 0.511945392491421};

	Matrix * A = new Matrix(m, n, data_A);
	Matrix * minusb = new Matrix(m, 1, data_minusb);
	Matrix * x0;
	Matrix xstar;
	Function * f = new QuadraticLoss();
	LinearOperator * OpA = new MatrixOperator(*A);
	Function * g = new Norm1(5.0);
	FBProblem prob = FBProblem(*f, *OpA, *minusb, *g);
	FBStoppingRelative sc = FBStoppingRelative(TOLERANCE);
	FBSplitting * solver;
	
	size_t repeat = 200;
	for (size_t r = 0; r < repeat; r++) {
		// test FB operations starting from x1
		x0 = new Matrix(n, 1, data_x1);
		solver = new FBSplittingFast(prob, *x0, gamma, sc, MAXIT);
		solver->run();
		xstar = solver->getSolution();
		//cout << "*** iters (fast) : " << solver->getIt() << endl;
		_ASSERT(solver->getIt() < MAXIT);
		for (int i=0; i < n; i++) {
			CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_xstar[i], xstar.get(i, 0), DOUBLES_EQUAL_DELTA);
		}
		delete x0;
		delete solver;
	}

	delete A;
	delete minusb;
	delete OpA;
	delete f;
	delete g;
}

void TestFBSplittingFast::testSparseLogReg_small() {
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
	double gamma = 0.1;
	// starting points
	double data_x1[] = {0, 0, 0, 0, 0};
	// reference results
	double ref_xstar[] = {0.0, 0.0, 0.215341883018748, 0.0, 0.675253988559914};

	Matrix * A = new Matrix(m, n, data_A);
	Matrix * minusb = new Matrix(m, 1, data_minusb);

	Function * f = new LogLogisticLoss();
	LinearOperator * OpA = new MatrixOperator(*A);
	Function * g = new Norm1(1.0);
	FBProblem prob(*f, *OpA, *minusb, *g);
	FBStoppingRelative sc(TOLERANCE);
	FBSplitting * solver;
	
	size_t repeat = 100;
	for (size_t r = 0; r < repeat; r++) {
		// test FB operations starting from x1
		Matrix * x0 = new Matrix(n, 1, data_x1);
		solver = new FBSplittingFast(prob, *x0, gamma, sc, MAXIT);
		solver->run();
		Matrix xstar = solver->getSolution();
		//cout << "*** iters (fast) : " << solver->getIt() << endl << flush;
		_ASSERT(solver->getIt() < MAXIT);
		for (int i=0; i < n; i++) {
			CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_xstar[i], xstar.get(i, 0), DOUBLES_EQUAL_DELTA);
		}
		delete x0;
		delete solver;
	}

	delete A;
	delete minusb;
	delete OpA;
	delete f;
	delete g;
}
