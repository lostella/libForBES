#ifndef ITERATIVESOLVER_H
#define	ITERATIVESOLVER_H

#include "Matrix.h"

/**
 * \class IterativeSolver
 * \brief Abstract class for iterative solvers
 */
class IterativeSolver {
protected:

    int m_it;
    int m_maxit;

    /**
     * Perform an iteration of the algorithm. As virtual method
     * its implementation is left to subclasses.
     *
     * \todo fix the return value.
     *
     * @return 0.
     */
    virtual int iterate() = 0;

    /**
     * Check stopping condition. As virtual method its implementation
     * is left to subclasses.
     *
     * @return 0 if the algorithm should NOT stop, nonzero otherwise.
     */
    virtual int stop() = 0;

public:

    /**
     * Constructs an iterative solver given the maximum number
     * of iterations to be executed.
     * 
     * @param maxit maximum number of iterations
     */
    explicit IterativeSolver(int maxit);

    /**
     * Runs the solver: until maximum number of iterations is met,
     * or stop() returns 1, execute iterate()
     *
     * \todo fix the return value
     */
    int run();

    /**
     * Gets the number of iterations performed by the solver.
     *
     * @return number of iterations performed.
     */
    int getIt();

    virtual ~IterativeSolver();

};

#endif /* ITERATIVESOLVER_H */