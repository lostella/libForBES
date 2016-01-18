#ifndef FBSPLITTING_H
#define	FBSPLITTING_H

#include "FBProblem.h"
#include "FBCache.h"
#include "IterativeSolver.h"

class FBSplitting : public IterativeSolver {
private:

    FBProblem * m_prob;
    double m_gamma;

    FBCache m_cache;
    double m_tol;

protected:

    virtual int iterate();
    virtual int stop();

public:

    /**
     * Initialize an FBSplitting object.
     *
     * @param p reference to the FBProblem to solve
     * @param x0 reference to Matrix, the starting point for the solver
     * @param gamma the initial stepsize parameter for the operations
     */
    FBSplitting(FBProblem & prob, Matrix & x0, double gamma);

    void setTol(double tol);

    /**
     * Gets the solution point computed by the algorithm.
     *
     * @return reference to a Matrix object, containing the problem solution
     */
    Matrix& getSolution();

    virtual ~FBSplitting();

};

#endif /* FBSPLITTING_H */
