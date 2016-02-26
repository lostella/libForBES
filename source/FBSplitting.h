#ifndef FBSPLITTING_H
#define	FBSPLITTING_H

#include "FBProblem.h"
#include "FBCache.h"
#include "IterativeSolver.h"
#include "FBStopping.h"

/**
 * \class FBSplitting
 * \brief Forward-backward splitting algorithm
 * \version version 0.0
 * \ingroup FBSolver-group
 * 
 * FBSplitting specializes IterativeSolver by defining how the
 * iterate() method works. In particular, it computes the next
 * iterate as the forward-backward (or proximal gradient) step
 * at the current point.
 */
class FBSplitting : public IterativeSolver {
private:

    /**
     * Specifications of the optimization problem.
     */
    FBProblem * m_prob;
    /**
     * Pointer to stopping criterion object.
     */
    FBStopping * m_sc;
    /**
     * Whether m_sc should be deleted in the class's destructor. This should be
     * set to \c true whenever a new FBStopping object is contructed internally
     * and to \c false whenever m_sc points to an external object.
     */
    bool delete_sc; 
    double m_gamma;

protected:

    FBCache m_cache;

public:

    virtual int iterate();
    virtual int stop();


public:

    /**
     * Initialize an FBSplitting object. By default, the maximum number
     * of iterations is set to 1000, and the tolerance on the fixed-point
     * residual is set to 1e-6.
     *
     * @param p reference to the FBProblem to solve
     * @param x0 reference to Matrix, the starting point for the solver
     * @param gamma the initial stepsize parameter for the operations
     */
    FBSplitting(FBProblem & prob, Matrix & x0, double gamma);

    /**
     * Initialize an FBSplitting object. By default, the maximum number
     * of iterations is set to 1000.
     *
     * @param p reference to the FBProblem to solve
     * @param x0 reference to Matrix, the starting point for the solver
     * @param gamma the initial stepsize parameter for the operations
     * @param sc reference to the FBStopping to be used as stopping criterion
     */
    FBSplitting(FBProblem & prob, Matrix & x0, double gamma, FBStopping & sc);

    /**
     * Initialize an FBSplitting object. By default, the tolerance on the
     * fixed-point residual is set to 1e-6.
     *
     * @param p reference to the FBProblem to solve
     * @param x0 reference to Matrix, the starting point for the solver
     * @param gamma the initial stepsize parameter for the operations
     * @param maxit maximum number of iterations
     */
    FBSplitting(FBProblem & prob, Matrix & x0, double gamma, int maxit);

    /**
     * Initialize an FBSplitting object.
     *
     * @param p reference to the FBProblem to solve
     * @param x0 reference to Matrix, the starting point for the solver
     * @param gamma the initial stepsize parameter for the operations
     * @param sc reference to the FBStopping to be used as stopping criterion
     * @param maxit maximum number of iterations
     */
    FBSplitting(FBProblem & prob, Matrix & x0, double gamma, FBStopping & sc, int maxit);

    /**
     * Gets the solution point computed by the algorithm.
     *
     * @return reference to a Matrix object, containing the problem solution
     */
    Matrix& getSolution();

    virtual ~FBSplitting();

};

#endif /* FBSPLITTING_H */
