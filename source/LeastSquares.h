/* 
 * File:   LeastSquares.h
 * Author: chung
 *
 * Created on January 14, 2016, 2:38 AM
 */

#ifndef LEASTSQUARES_H
#define	LEASTSQUARES_H

#include "MatrixSolver.h"
#include "lapacke.h"


class LeastSquares : public MatrixSolver {
public:
    LeastSquares(Matrix& matrix);
    
    virtual ~LeastSquares();
    
    virtual int solve(Matrix& rhs, Matrix& solution);
    
    virtual int solve(Matrix& rhs_in_out);

private:

};

#endif	/* LEASTSQUARES_H */

