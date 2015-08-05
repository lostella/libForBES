/* 
 * File:   LDLFactorization.h
 * Author: chung
 *
 * Created on July 30, 2015, 3:02 AM
 */

#ifndef LDLFACTORIZATION_H
#define	LDLFACTORIZATION_H

#include <cstring>

#include "Matrix.h"
#include "FactoredSolver.h"

class LDLFactorization : public FactoredSolver {
public:

    LDLFactorization(Matrix& m_matrix);

    virtual ~LDLFactorization();

    virtual int factorize(void);
    
    virtual int solve(const Matrix& rhs, Matrix& solution) const;
    
    double* getLDL() const;

    int* getIpiv() const;

private:


    double* LDL = NULL;
    int* ipiv = NULL;

};

#endif	/* LDLFACTORIZATION_H */

