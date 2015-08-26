/* 
 * File:   FactoredSolver.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 30, 2015, 3:18 AM
 */

#ifndef FACTOREDSOLVER_H
#define	FACTOREDSOLVER_H

#include "Matrix.h"

class FactoredSolver {
public:

    FactoredSolver(Matrix& m_matrix) : m_matrix(m_matrix) {}

    virtual ~FactoredSolver();

    virtual int factorize(void) = 0;
    virtual int solve(const Matrix& rhs, Matrix& solution) const = 0;

private:

protected:

    Matrix& m_matrix;

};

#endif	/* FACTOREDSOLVER_H */

