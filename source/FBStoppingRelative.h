#ifndef FBSTOPPINGRELATIVE_H
#define	FBSTOPPINGRELATIVE_H

#include "FBStopping.h"
#include "FBCache.h"

/**
 * @brief Relative stopping criterion on the fixed-point residual
 * 
 * This class overrides the stop(FBCache& c) method of FBStopping
 * by checking the Euclidean norm of the fixed-point residual relative
 * to the norm of the current iterate. That is, given
 * \f$\mbox{tol}\f$, the algorithm will stop if
 * 
 * \f[ \|x - T_\gamma(x)\|_2 \leq \mbox{tol}(1+\|x\|_2), \f]
 * 
 * where \f$x\f$ is the current iterate.
 */
class FBStoppingRelative : public FBStopping {

public:

    FBStoppingRelative(double tol);

	virtual int stop(FBCache & c);

	virtual ~FBStoppingRelative();
	
};

#endif /* FBSTOPPINGRELATIVE_H */
