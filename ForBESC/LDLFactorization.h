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

    LDLFactorization(Matrix& m_matrix) :
    FactoredSolver(m_matrix) {
        if (m_matrix.getType() != Matrix::MATRIX_DENSE) {
            throw std::logic_error("This matrix type is not supported by LDLFactorization");
        }
        if (m_matrix.getNrows() != m_matrix.getNcols()) {
            throw std::logic_error("Matrix not square");
        }
        LDL = new double[m_matrix.length()];
        ipiv = new int[m_matrix.getNrows()];
        memcpy(LDL, m_matrix.getData(), m_matrix.length() * sizeof (double));
    }


    virtual ~LDLFactorization();

    int factorize(void) {
        size_t n = m_matrix.getNrows();
        int status = LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'L', n, LDL, n, ipiv);
        return status;
    }

    virtual int solve(const Matrix& rhs, Matrix& solution) {
        size_t n = m_matrix.getNrows();
        if (solution.length() < n) {
            solution = Matrix(rhs);
        } else if (solution.length() < n) {
            solution.reshape(n, 1);
        }
        int status = LAPACKE_dsytrs(LAPACK_COL_MAJOR, 'L', n, 1, LDL, n, ipiv, solution.getData(), n);
        return status;
    }

    double* getLDL() const {
        return LDL;
    }

    int* getIpiv() const {
        return ipiv;
    }

private:


    double* LDL = NULL;
    int* ipiv = NULL;

};

#endif	/* LDLFACTORIZATION_H */

