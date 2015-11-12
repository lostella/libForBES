#ifndef SOLVER_H
#define	SOLVER_H

class Solver {
public:

	int run() {
		while(!stop()) {
			prepare();
			iterate();
		}
	}

	virtual ~Solver();

protected:

	virtual int prepare() = 0;
	virtual int iterate() = 0;
	virtual int stop() = 0;

};

#include "FBProblem.h"

class FBSolver : public Solver {
public:

	int setMaxIters(int maxIt);
	int setTolerance(double tol);

protected:

  FBProblem prob;

  FBSolver(FBProblem & prob) : prob(prob) {}

  virtual ~FBSolver() {}

private:

	int it, maxIt;
	double tol;

};

#endif /* SOLVER_H */
