/* 
 * File:   LDLFactorization_AAt.h
 * Author: chung
 *
 * Created on November 5, 2015, 1:12 AM
 */

#ifndef LDLFACTORIZATION_AAT_H
#define	LDLFACTORIZATION_AAT_H

#include "FactoredSolver.h"

/**
 * 
 * \brief Cholesky factorization and solver.
 */
class S_LDLFactorization : public FactoredSolver {
public:
    
    S_LDLFactorization(Matrix& matrix, double beta);

    virtual ~S_LDLFactorization();

    virtual int factorize();

    virtual int solve(const Matrix& rhs, Matrix& solution) const;

private:

    cholmod_factor * m_factor;
    double m_beta;

};

#endif	/* LDLFACTORIZATION_AAT_H */

