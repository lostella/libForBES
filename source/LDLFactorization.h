/* 
 * File:   LDLFactorization.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 30, 2015, 3:02 AM
 */

#ifndef LDLFACTORIZATION_H
#define	LDLFACTORIZATION_H

#include <cstring>

#include "ForBESUtils.h"
#include "Matrix.h"
#include "FactoredSolver.h"
#include "ldl.h"

/**
 * 
 * \brief LDL factorization and solver.
 */
class LDLFactorization : public FactoredSolver {
public:

    LDLFactorization(Matrix& m_matrix);

    virtual ~LDLFactorization();

    virtual int factorize(void);

    virtual int solve(const Matrix& rhs, Matrix& solution) const;

    double* getLDL() const;

    int* getIpiv() const;

private:


    double* LDL = NULL; /* LDL factorization computed by lapack */
    int* ipiv = NULL; /* Pivots for the LDL factorization computed by lapack */

    /**
     * A sparse LDL factorization
     */
    typedef struct sparse_ldl_factor_struct {
        double * Lx; /**< Values of L */
        int * Li; /**< i-pointers of L */
        int * Lp;/**< p-pointers of L */
        double * D; /**< Diagonal part of the LDL factorization*/
    } sparse_ldl_factor;

    /**
     * Pointer to a sparse LDL factorization.
     */
    sparse_ldl_factor * m_sparse_ldl_factor = NULL;

    

};

#endif	/* LDLFACTORIZATION_H */

