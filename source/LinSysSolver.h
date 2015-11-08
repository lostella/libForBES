/* 
 * File:   LinSysSolver.h
 * Author: chung
 *
 * Created on November 7, 2015, 3:47 PM
 */

#ifndef LINSYSSOLVER_H
#define	LINSYSSOLVER_H

#include "Matrix.h"
/**
 * \class LinSysSolver
 * \brief An abstract solver for linear systems
 * \version version 0.0
 * \date Created on November 7, 2015, 3:47 PM
 * \author Pantelis Sopasakis
 * 
 * LinSysSolver is a solver interfact for linear systems of the form \f$Ax=b\f$.
 */
class LinSysSolver {
public:
        
    
    virtual ~LinSysSolver();

    /**
     * Solves the linear system.
     * 
     * @param rhs the right-hand side of the linear equation
     * @param solution the solution of the linear system which is computed using the 
     * matrix factorization.
     * @return status code. The method returns \link ForBESUtils::STATUS_OK STATUS_OK\endlink 
     * if the invocation has succeeded or 
     * \link ForBESUtils::STATUS_NUMERICAL_PROBLEMS STATUS_NUMERICAL_PROBLEMS\endlink
     * if numerical problems have hindered the computation of a solution.
     * 
     */
    virtual int solve(Matrix& rhs, Matrix& solution) const = 0;

private:

protected:
    
    /**
     * Instantiate a new LinSysSolver object
     * @param matrix
     */
    explicit LinSysSolver(Matrix& matrix);

    Matrix& m_matrix;
    Matrix::MatrixType m_matrix_type;
    size_t m_matrix_nrows;
    size_t m_matrix_ncols;

};

#endif	/* LINSYSSOLVER_H */

