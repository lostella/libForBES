#ifndef FBSPLITTING_H
#define	FBSPLITTING_H

#include "FBProblem.h"
#include "Solver.h"
#include "Solution.h"

class FBSplitting : public FBSolver {
private:

	Matrix x, z;

public:

	FBSplitting(FBProblem & prob, Matrix & x0, double gamma) :
		FBSolver(prob), this->prob(prob), this->x(x0), this.gamma(gamma) {}

	virtual ~FBSplitting();

};

#endif /* FBSPLITTING_H */
