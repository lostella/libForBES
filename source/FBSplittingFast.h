#ifndef FBSPLITTINGFAST_H
#define	FBSPLITTINGFAST_H

#include "FBProblem.h"
#include "FBCache.h"
#include "IterativeSolver.h"
#include "FBSplitting.h"
#include "FBStopping.h"

/**
 * \class FBSplittingFast
 * \brief Fast forward-backward splitting algorithm
 * \version version 0.0
 * \ingroup FBSolver-group
 * 
 * FBSplittingFast extends FBSplitting by prepending the Nesterov
 * extrapolation step to the <code>iterate()</code> method. On convex problems
 * this is known to enforce the optimal convergence rate (for the
 * objective value) of order \f$\mathcal{O}(1/k^2)$\f$.
 */
class FBSplittingFast : public FBSplitting {
private:

    Matrix * m_previous;

protected:

public:

    virtual int iterate();
    virtual int stop();

    /**
     * Initialize an FBSplittingFast object. By default, the maximum number
     * of iterations is set to \c 1000, and the tolerance on the fixed-point
     * residual is set to \c 1e-6.
     *
     * @param p reference to the FBProblem to solve
     * @param x0 reference to Matrix, the starting point for the solver
     * @param gamma the initial stepsize parameter for the operations
     */
    FBSplittingFast(FBProblem & prob, Matrix & x0, double gamma);

    /**
     * Initialize an FBSplittingFast object. By default, the maximum number
     * of iterations is set to \c 1000.
     *
     * @param p reference to the FBProblem to solve
     * @param x0 reference to Matrix, the starting point for the solver
     * @param gamma the initial stepsize parameter for the operations
     * @param sc reference to the FBStopping to be used as stopping criterion
     */
    FBSplittingFast(FBProblem & prob, Matrix & x0, double gamma, FBStopping & sc);

    /**
     * Initialize an FBSplittingFast object. By default, the tolerance on the
     * fixed-point residual is set to \c 1e-6.
     *
     * @param p reference to the FBProblem to solve
     * @param x0 reference to Matrix, the starting point for the solver
     * @param gamma the initial stepsize parameter for the operations
     * @param maxit maximum number of iterations
     */
    FBSplittingFast(FBProblem & prob, Matrix & x0, double gamma, int maxit);

    /**
     * Initialize an FBSplittingFast object.
     *
     * @param p reference to the FBProblem to solve
     * @param x0 reference to Matrix, the starting point for the solver
     * @param gamma the initial stepsize parameter for the operations
     * @param sc reference to the FBStopping to be used as stopping criterion
     * @param maxit maximum number of iterations
     */
    FBSplittingFast(FBProblem & prob, Matrix & x0, double gamma, FBStopping & sc, int maxit);

    /**
     * Gets the solution point computed by the algorithm.
     *
     * @return reference to a Matrix object, containing the problem solution
     */
    Matrix& getSolution();

    virtual ~FBSplittingFast();

};

#endif /* FBSPLITTINGFAST_H */
