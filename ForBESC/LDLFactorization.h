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

class LDLFactorization {
public:

    LDLFactorization(Matrix& matrix) :
    matrix(matrix) {
        if (matrix.getType() != Matrix::MATRIX_DENSE) {
            throw std::logic_error("This matrix type is not supported by LDLFactorization");
        }
        if (matrix.getNrows() != matrix.getNcols()) {
            throw std::logic_error("Matrix not square");
        }
        LDL = new double[matrix.length()];
        ipiv = new int[matrix.getNrows()];
        memcpy(LDL, matrix.getData(), matrix.length() * sizeof (double));
    }

    virtual ~LDLFactorization();

    int factorize(void) {
        size_t n = matrix.getNrows();
        int status = LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'L', n, LDL, n, ipiv);
        return status;
    }

    double* getLDL() const {
        return LDL;
    }

    int* getIpiv() const {
        return ipiv;
    }

private:

    Matrix& matrix;
    double* LDL = NULL;
    int* ipiv = NULL;

};

#endif	/* LDLFACTORIZATION_H */

