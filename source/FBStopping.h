#ifndef FBSTOPPING_H
#define	FBSTOPPING_H

#include "FBCache.h"

/**
 * @brief Basic stopping criterion for forward-backward problems
 * 
 * This class implements the most basic stopping criterion based
 * on the norm of the fixed-point residual. Given the tolerance
 * \f$\mbox{tol}\f$, the algorithm will stop if
 * 
 * \f[ \|x - T_\gamma(x)\|_2 \leq \mbox{tol}, \f]
 * 
 * where \f$x\f$ is the current iterate.
 */
class FBStopping {
protected:

    double m_tol;

public:

    /**
     * @brief Instantiate FBStopping given a tolerance
     * 
     * @param tol the prescribed tolerance on the fixed-point
     * residual.
     */
    FBStopping(double tol);

    /**
     * @brief Checks whether the algorithm should stop
     * 
     * @param c FBCache object containing the current state of
     * the algorithm.
     * @return 0 if the algorithm should go on, nonzero if it
     * should stop.
     */
    virtual int stop(FBCache & c);

    virtual ~FBStopping();

};

#endif /* FBSTOPPING_H */
