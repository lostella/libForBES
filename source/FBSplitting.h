#ifndef FBSPLITTING_H
#define	FBSPLITTING_H

#include "FBProblem.h"
#include "Solver.h"
//#include "Solution.h"

class FBSplitting : public FBSolver {
private:

    Matrix x, z;
    double gamma;
public:

    FBSplitting(FBProblem & prob, Matrix & x0, double gamma) :
    FBSolver(prob),
    x(x0),
    gamma(gamma) {
    }

    virtual ~FBSplitting();

};

#endif /* FBSPLITTING_H */
