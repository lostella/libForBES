/* 
 * File:   Matrix.h
 * Author: Chung
 *
 * Created on July 7, 2015, 8:02 PM
 * 
 * ForBES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *  
 * ForBES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with ForBES. If not, see <http://www.gnu.org/licenses/>.
 */

//TODO Implement packed storage for symmetric and lower/upper triangular matrices
//TODO Implement multiplication with such matrices
//TODO Use CSparse to introduce sparse matrices

#ifndef MATRIX_H
#define	MATRIX_H

#include <iostream>
#include <stdexcept>
#include <complex>
#include <cstdlib>
#include <stdlib.h>
#include <iomanip>
#include <cmath>

#include "cholmod.h"

#ifdef USE_LIBS
#include <cblas.h>
#include <lapacke.h>
#endif

/**
 * A m-by-n matrix (or a vector).
 */
class Matrix {
public:

    /**
     * Types of matrices.
     */
    enum MatrixType {
        MATRIX_DENSE,
        MATRIX_SPARSE,
        MATRIX_DIAGONAL,
        MATRIX_LOWERTR,
        MATRIX_SYMMETRIC
    };

    /**
     * Sparse type of matrix
     */
    enum SparseMatrixType {
        SPARSE_UNSYMMETRIC = 0,
        SPARSE_SYMMETRIC_L = 1,
        SPARSE_SYMMETRIC_R = -1
    };

    /* Constructors and destructors */

    /**
     * Default constructor for <code>Matrix</code>.
     */
    Matrix();

    /**
     * 
     * @param nr
     * @param nc
     */
    Matrix(int nr, int nc);

    /**
     * 
     * @param nr
     * @param nc
     * @param matrixType
     */
    Matrix(int nr, int nc, MatrixType matrixType);

    /**
     * 
     * @param nr
     * @param nc
     * @param data
     */
    Matrix(int nr, int nc, const float * data);

    /**
     * 
     * @param nr
     * @param nc
     * @param data
     * @param matrixType
     */
    Matrix(int nr, int nc, const float * data, MatrixType matrixType);

    /**
     * Copy-constructor.
     * @param orig
     */
    Matrix(const Matrix& orig);

    /**
     * Destructor.
     * Sets the row and column dimension of this matrix to <code>-1</code> and
     * deletes the matrix data.
     */
    virtual ~Matrix();



    /**
     * Transpose this matrix. 
     */
    void transpose();


    /**
     * Returns the value of the current matrix at position <code>(i,j)</code>.
     * @param i row index (<code>0,...,nrows-1</code>)
     * @param j column index (<code>0,...,ncols-1</code>)
     * @return value at <code>(i,j)</code>
     */
    float get(const int i, const int j) const;


    /**
     * Updates a value of the matrix at position <code>(i,j)</code>
     * @param i row index (<code>0,...,nrows-1</code>)
     * @param j column index (<code>0,...,ncols-1</code>)
     * @param val value to be set at <code>(i,j)</code>
     */
    void set(int i, int j, float val);


    /* Getters */
    /**
     * Get the number of columns of the current matrix.
     * @return columns as <code>int</code>.
     */
    int getNcols() const;

    /**
     * Get the number of rows of the current matrix.
     * @return rows as <code>int</code>.
     */
    int getNrows() const;

    /**
     * Getter for the matrix data. Provides direct access to the matrix data which
     * are stored as <code>float *</code>.
     * 
     * @return Pointer to the matrix data
     */
    float * const getData() const;

    /**
     * Returns the type of this matrix as <code>MatrixType</code>
     * @return 
     */
    MatrixType getType() const;

    /**
     * Reshape the matrix.
     * @param nrows new number of rows
     * @param ncols new number of columns
     * @return status code: <code>0</code> if reshaping succeeded, <code>-1</code>
     * if some of the new dimensions is 0, <code>-2</code> if reshaping is 
     * impossible.
     */
    int reshape(int nrows, int ncols);



    /* Utilities */

    /**
     * 
     * @return <code>true</code> if this is a column-vector and <code>false</code>
     * otherwise. 
     */
    bool isColumnVector() const;

    /**
     * 
     * @return <code>true</code> if this is a row-vector and <code>false</code>
     * otherwise.
     */
    bool isRowVector() const;

    /**
     * 
     * @return <code>true</code> if this is an empty vector.
     */
    bool isEmpty() const;

    /**
     * Length of data of this matrix (e.g., if this is a diagonal matrix, only its
     * diagonal elements are stored, so the data length equals the row-dimension 
     * of the matrix).
     * @return Data length.
     */
    int length() const;

    /**
     * Computes the quadratic form x'*Q*x, where x is a given vector
     * as <code>Matrix</code>, where Q is the current instance of 
     * <code>Matrix</code>.
     * 
     * This method can only be applied on square matrices Q while x 
     * needs to be of appropriate dimension.
     * 
     * @param x The vector x.
     * @return Scalar x'*Q*x as <code>float</code>.
     */
    float quad(const Matrix& x) const;

    /**
     * Computes the quadratic form x'*Q*x + q'*x, where x is a given vector
     * as <code>Matrix</code>, where Q is the current instance of 
     * <code>Matrix</code>. Computes x'*Q*x.
     * 
     * This method can only be applied on square matrices Q while x and q 
     * need to be of appropriate dimensions.
     * 
     * @param x The vector x.
     * @param q The parameter vector x.
     * @return 
     */
    float quad(const Matrix& x, const Matrix& q) const;

    /**
     * Computes the Cholesky factorization of this matrix. If this is a <code>DENSE<code>
     * matrix, then it is assumed it is symmetric (but there is no verification) and
     * only its lower triangular part is considered. Notice that the Cholesky factorization
     * can only be applied to symmetric and poisitive definite matrices.
     * 
     * This method, when applied on a <code>DENSE<code> matrix, the produced matrix <code>L</code>
     * will be of type <code>DENSE<code>. It is advisable to apply this method only
     * on matrices of type <code>SYMMETRIC</code>.
     * 
     * @param L the cholesky factor of this matrix. 
     * @return status code. Returns <code>0</code> if the factorization succeeded.
     */
    int cholesky(Matrix& L);

    /**
     * Solves the linear system <code>L*L'*x = b</code>, where <code>L</code> is the
     * current matrix, <code>b</code> is a given right-hand side vector or matrix
     * (given as an instance of <code>Martrix</code>) and <code>x</code> is the
     * solution which is returned by this method.
     * 
     * Note that if this is a <code>DENSE<code> matrix, then it is assumed it is 
     * lower triangular (but there is no verification) and only its lower triangular 
     * part is considered.
     * 
     * @param solution the solution of the linear system as an instance of <code>Matrix</code>
     * @param rhs The right-hand side vector or matrix
     * 
     * @return Returns <code>0</code> if the solution of the system has succeeded.
     */
    int solveCholeskySystem(Matrix& solution, const Matrix& rhs) const;


    //TODO Implement LDL factorization
    //TODO Implement LDL-based system solution

    /* Operators */

    /**
     * 
     * @param sub index
     * @return 
     */
    float &operator[](const int sub) const; //overloading []    

    /**
     * 
     * @param right is the right-hand side matrix
     * @return 
     */
    Matrix operator+(const Matrix& right) const;

    /**
     * 
     * @param right is the right-hand side matrix
     * @return 
     */
    Matrix operator-(const Matrix& right) const;

    /**
     * 
     * @param right is the right-hand side matrix
     * @return 
     */
    Matrix& operator+=(const Matrix& right);

    /**
     * 
     * @param right is the right-hand side matrix
     * @return 
     */
    Matrix& operator-=(const Matrix& right);

    /**
     * Overloaded multiplication operator for <code>Matrix</code>
     * @param right is the right-hand side matrix
     * @return 
     */
    Matrix operator*(const Matrix& right) const;

    /**
     * Assignment operator.
     * @param right is the right-hand operand.
     * @return 
     */
    Matrix& operator=(const Matrix& right);

    /**
     * Equality relational operator: returns <code>true</code> iff both sides
     * are equal. Two matrices are equal if they are of the same type, have equal
     * dimensions and equal values.
     * @param right
     * @return 
     */
    bool operator==(const Matrix& right) const;


    /**
     * @param os
     * @param obj
     */
    friend std::ostream& operator<<(std::ostream& os, const Matrix& obj);

    cholmod_common  *m_cholmod_common = NULL;
    cholmod_triplet *m_triplet = NULL;
    cholmod_sparse  *m_sparse = NULL;
    cholmod_factor  *m_cholesky_factor = NULL;
    
private:
    /*
     * MatrixFactory is allowed to access these private fields!
     */
    friend class MatrixFactory;
    
    int m_nrows; /*< Number of rows */
    int m_ncols; /*< Number of columns */
    int m_dataLength; /*< Length of data */
    float *m_data = NULL; /*< Data */
    bool m_transpose;
    MatrixType m_type; /*< Matrix type */
    
    
    void createSparseFromTriplet();

    /**
     * Initialize the current matrix (allocate memory etc) for a given number of 
     * rows and columns and a given matrix type.
     * @param nrows Number of rows
     * @param ncols Number of column
     * @param matrixType Matrix type
     */
    void init(int nrows, int ncols, MatrixType matrixType);

    /**
     * Check whether a given pair of indexes is within the matrix bounds.
     * @param i row index
     * @param j column 
     * @return <code>true</code> if (i,j) is within the bounds
     */
    bool indexWithinBounds(int i, int j);

    /**
     * Multiply with a matrix when the left-hand side matrix is dense
     * @param right any right-hand side matrix
     * @return 
     */
    Matrix multiplyLeftDense(const Matrix& right) const;

    /**
     * Multiply with a matrix when the left-hand side matrix is diagonal.
     * @param right any right-hand side matrix
     * @return 
     */
    Matrix multiplyLeftDiagonal(const Matrix& right) const;

    Matrix multiplyLeftSymmetric(const Matrix& right) const;

    void domm(const Matrix &right, Matrix &result) const;

    enum SparseMatrixStorageType {
        CHOLMOD_TYPE_TRIPLET = 444,
        CHOLMOD_TYPE_SPARSE = 555,
        CHOLMOD_TYPE_DENSE = 666,
        CHOLMOD_TYPE_FACTOR = 777
    };
    
    SparseMatrixStorageType m_sparseStorageType;
};

#endif	/* MATRIX_H */

