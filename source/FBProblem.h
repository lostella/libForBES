#ifndef FBPROBLEM_H
#define	FBPROBLEM_H

#include "Function.h"
#include "LinearOperator.h"

class FBProblem {
private:
  Function * f1, * f2, * g;
  LinearOperator * L1, * L2;
  Matrix * d1, * d2, * lin;

public:
  FBProblem() {}
	FBProblem(Function * f1, LinearOperator * L1, Matrix * d1,
    Function * f2, LinearOperator * L2, Matrix * d2, Matrix * lin,
    Function * g);

	Function        * getf1() { return f1; }
	LinearOperator  * getL1() { return L1; }
	Matrix          * getd1() { return d1; }
	Function        * getf2() { return f2; }
	LinearOperator  * getL2() { return L2; }
	Matrix          * getd2() { return d2; }
	Matrix          * getlin() { return lin; }
	Function        * getg()  { return g; }

	virtual ~FBProblem() {}

};

#endif /* FBPROBLEM_H */
