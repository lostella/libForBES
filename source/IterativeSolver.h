#ifndef ITERATIVESOLVER_H
#define	ITERATIVESOLVER_H

#include "Matrix.h"

class IterativeSolver {
protected:

	virtual int iterate() = 0;
	virtual int stop() = 0;

public:

	int run();

	virtual ~IterativeSolver();
	
};

#endif /* ITERATIVESOLVER_H */

