#include "FBProblem.h"

#include "TestLasso.h"

#define DOUBLES_EQUAL_DELTA 1e-8
#define MAXIT 1000
#define TOLERANCE 1e-12

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
	/* Define the problem data */
	const size_t n = 5;
	const size_t m = 4;
	double data_A[] = {
		1,  2, -1, -1,
		-2, -1,  0, -1,
		3,  0,  4, -1,
		-4, -1, -3,  1,
		5,  3,  2,  3 };
	double data_minus_b[] = {-1, -2, -3, -4};
	Matrix A(m, n, data_A);
	Matrix minus_b(m, 1, data_minus_b);

	double ref_xstar[] = {-0.010238907850120, 0.0, 0.0, 0.0, 0.511945392491510};

	LinearOperator * OpA = new MatrixOperator(A);
	Function * f = new QuadraticLoss();
	double lambda = 5.0;
	Function * g = new Norm1(lambda);

	FBStoppingRelative sc(TOLERANCE);

	// Define the FB problem
	FBProblem prob = FBProblem(*f, *OpA, minus_b, *g);
	// Initial guess and gamma - Construct a new instance of FBSplitting
	Matrix x0(n, 1);
	double gamma = 1e-2;
	FBSplitting * solver = new FBSplitting(prob, x0, gamma, sc, MAXIT);
	// Run the solver and get the solution
	solver->run();
	Matrix xstar = solver->getSolution();

	for (int i=0; i < n; i++) {
		CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_xstar[i], xstar.get(i, 0), DOUBLES_EQUAL_DELTA);
	}
}