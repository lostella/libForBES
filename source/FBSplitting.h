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
    int m_it;
    int m_maxit;
    double m_tol;

protected:

	virtual int iterate();
	virtual int stop();
    
public:

    FBSplitting(FBProblem & prob, Matrix & x0, double gamma);
    
    int setMaxIt(int maxit);
    int setTol(double tol);
	
	int getIt();
	Matrix& getSolution();

    virtual ~FBSplitting();

};

#endif /* FBSPLITTING_H */
