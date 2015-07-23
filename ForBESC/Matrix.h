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

//TODO Implement packed storage for symmetric and lower/upper triangular matrices [done]
//TODO Implement multiplication with such matrices [almost done]
//TODO Use CSparse to introduce sparse matrices [ongoing]

#ifndef MATRIX_H
#define	MATRIX_H

#include <iostream>
#include <stdexcept>
#include <complex>
#include <cstdlib>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <cstring>

#include "cholmod.h"

#ifdef USE_LIBS
#include <cblas.h>
#include <lapacke.h>
#endif

/**
 * A generic matrix API.
 * 
 * <p>A generic matrix which can be an unstructured dense, a structured dense (e.g.,
 * symmetric or lower triangular, stored in packed form) or a dense matrix. 
 * This class provides a uniform access framework (an API) to matrix-matrix
 * operations (e.g., addition and multiplication), factorizations and other useful 
 * operations.</p>
 * 
 * <p>To construct a Matrix you can use one of this class's constructors. However,
 * for sparse matrices it is advisable to use the factory class <code>MatrixFactory</code>.</p>
 */
class Matrix {
    
public:
    
    /* STATIC */
    
    /**
     * This is the single access method to the singleton <code>cholmod_common</code>
     * used in this project. Typically clients will not be interested in using this
     * <code>cholmod_common</code> to perform any matrix-matrix operations or factorization,
     * however, it can be used to check the status of computations, get the overall
     * flop count and more. 
     * 
     * This method will construct and store internally an instance of <code>cholmod_common</code>
     * if one does not exist.
     * 
     * @return The singleton <code>cholmod_common</code> object.
     */
    static cholmod_common* cholmod_handle();
    
    static int destroy_handle();

    /**
     * Types of matrices.
     */
    enum MatrixType {
        MATRIX_DENSE,       /**< A dense matrix */
        MATRIX_SPARSE,      /**< A sparse matrix (powered by SuiteSparse) */
        MATRIX_DIAGONAL,    /**< A diagonal matrix */
        MATRIX_LOWERTR,     /**< A lower-triangular matrix */
        MATRIX_SYMMETRIC    /**< A symmetric matrix */
    };

    /**
     * Sparse type of matrix
     */
    enum SparseMatrixType {
        SPARSE_UNSYMMETRIC = 0, /**< Not symmetric sparse */
        SPARSE_SYMMETRIC_L = 1, /**< Symmetric sparse (lower-triangular part is used) */
        SPARSE_SYMMETRIC_R = -1 /**< Symmetric sparse (upper-triangular part is used) */
    };

    /* Constructors and destructors */

    /**
     * Default constructor for <code>Matrix</code>.
     */
    Matrix();

    /**
     * Allocates an empty dense matrix of given dimensions
     * 
     * @param nr number of rows
     * @param nc number of columns
     */
    Matrix(size_t nr, size_t nc);

    /**
     * Allocates a matrix of given dimensions and given type.
     * 
     * @param nr number of rows
     * @param nc number of columns
     * @param matrixType type of matrix
     */
    Matrix(size_t nr, size_t nc, MatrixType matrixType);

    /**
     * Allocates a new dense matrix.
     * 
     * @param nr number of rows
     * @param nc number of columns
     * @param data double values (data will be copied)
     */
    Matrix(size_t nr, size_t nc, const double * data);

    /**
     * Allocates a new matrix of given dimensions, given data and given 
     * matrix type. Use this constructor only for non-sparse matrices; it is
     * recommended to use the factory class <code>MatrixFactory</code> to 
     * construct instances of sparse matrices.
     * 
     * @param nr number of rows
     * @param nc number of columns
     * @param data double values (data will be copied)
     * @param matrixType a non-sparse matrix type.
     */
    Matrix(size_t nr, size_t nc, const double * data, MatrixType matrixType);

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
    double get(const size_t i, const size_t j) const;


    /**
     * Updates a value of the matrix at position <code>(i,j)</code>
     * @param i row index (<code>0,...,nrows-1</code>)
     * @param j column index (<code>0,...,ncols-1</code>)
     * @param val value to be set at <code>(i,j)</code>
     */
    void set(size_t i, size_t j, double val);


    /* Getters */

    /**
     * Get the number of columns of the current matrix.
     * @return columns as <code>int</code>.
     */
    size_t getNcols() const;

    /**
     * Get the number of rows of the current matrix.
     * @return rows as <code>int</code>.
     */
    size_t getNrows() const;

    /**
     * Getter for the matrix data. Provides direct access to the matrix data which
     * are stored as an array of <code>double</code> (datatype <code>double*</code>).
     * 
     * @return Pointer to the matrix data
     */
    double * const getData() const;

    /**
     * Returns the type of this matrix as <code>MatrixType</code>
     * @return 
     */
    MatrixType getType() const;

    /**
     * Reshape the matrix.
     * 
     * @param nrows new number of rows
     * @param ncols new number of columns
     * 
     * @return status code: <code>0</code> if reshaping succeeded, <code>-1</code>
     * if some of the new dimensions is 0, <code>-2</code> if reshaping is 
     * impossible.
     */
    int reshape(size_t nrows, size_t ncols);



    /* Utilities */

    /**
     * Checks whether this is a column vector. 
     * 
     * @return <code>true</code> if this is a column-vector and <code>false</code>
     * otherwise. 
     */
    bool isColumnVector() const;

    /**
     * Checks whether this is a row vector.
     * 
     * @return <code>true</code> if this is a row-vector and <code>false</code>
     * otherwise.
     */
    bool isRowVector() const;

    /**
     * Checks whether this is an empty matrix.
     * 
     * @return <code>true</code> if this is an empty vector.
     */
    bool isEmpty() const;

    /**
     * Length of data of this matrix (e.g., if this is a diagonal matrix, only its
     * diagonal elements are stored, so the data length equals the row-dimension 
     * of the matrix).
     * 
     * @return Data length.
     */
    size_t length() const;

    /**
     * Computes the quadratic form x'*Q*x.
     * 
     * <p>Here x is a given vector
     * as <code>Matrix</code>, where Q is the current instance of 
     * <code>Matrix</code>.</p>
     * 
     * <p>This method can only be applied on square matrices Q while x and q 
     * need to be of compatible dimensions.</p>
     * 
     * @param x The vector x.
     * @return Scalar x'*Q*x as <code>double</code>.
     */
    double quad(Matrix& x);

    /**
     * Computes the quadratic form x'*Q*x + q'*x.
     * 
     * <p>Here x is a given vector
     * as <code>Matrix</code>, where Q is the current instance of 
     * <code>Matrix</code>. Computes x'*Q*x.</p>
     * 
     * <p>This method can only be applied on square matrices Q while x and q 
     * need to be of combatible dimensions.</p>
     * 
     * @param x The vector x.
     * @param q The parameter vector x.
     * @return The result of x'*Q*x + q'*x.
     */
    double quad(Matrix& x, Matrix& q);

    /**
     * Computes the Cholesky factorization of this matrix. 
     * 
     * <p>If this is a <code>MATRIX_DENSE</code>
     * matrix, then it is assumed it is symmetric (but there is no verification) and
     * only its lower triangular part is considered. Notice that the Cholesky factorization
     * can only be applied to symmetric and positive definite matrices.</p>
     * 
     * <p>This method, when applied on a <code>MATRIX_DENSE</code> matrix, the produced matrix <code>L</code>
     * will be of type <code>MATRIX_DENSE</code>. It is advisable to apply this method only
     * on matrices of type <code>MATRIX_SYMMETRIC</code>.</p>
     * 
     * <p>Note that this is not a <code>const</code> method since, especially in the
     * case of sparse matrices, it may update the internal state of the matrix to 
     * facilitate and speed-up computations</p>
     * 
     * @param L the Cholesky factor of this matrix. 
     * @return status code. Returns <code>0</code> if the factorization succeeded.
     * 
     */
    int cholesky(Matrix& L);

    /**
     * Solves the linear system <code>L*L'*x = b</code>, where <code>L</code> is the
     * current matrix, <code>b</code> is a given right-hand side vector or matrix
     * (given as an instance of <code>Martrix</code>) and <code>x</code> is the
     * solution which is returned by this method.
     * 
     * Note that if this is a <code>MATRIX_DENSE</code> matrix, then it is assumed it is 
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
     * Direct access to the matrix data.
     * Shorthand for <code>matrix.getData()[]</code>; it is however safer to access
     * the matrix entries using <code>get</code> and <code>set</code>.
     * 
     * @param sub index
     * @return reference to matrix data
     */
    double &operator[](const int sub) const; //overloading []    

    /**
     * Summation operator.
     * 
     * Matrices must have compatible dimensions. 
     * @param right is the right-hand side matrix
     * @return 
     */
    Matrix operator+(Matrix& right) const;

    /**
     * Subtraction operator. 
     * 
     * Matrices must have compatible dimensions. 
     * @param right is the right-hand side matrix
     * @return 
     */
    Matrix operator-(const Matrix& right) const;

    /**
     * 
     * @param right is the right-hand side matrix
     * @return updated instance of <code>Matrix</code> 
     */
    Matrix& operator+=(Matrix& right);

    /**
     * 
     * @param right is the right-hand side matrix
     * @return updated instance of <code>Matrix</code>
     */
    Matrix& operator-=(const Matrix& right);

    /**
     * Overloaded multiplication operator for <code>Matrix</code>.
     * @param right is the right-hand side matrix
     * @return the result of the multiplication of two matrices
     */
    Matrix operator*(Matrix& right);

    /**
     * Assignment operator.
     * @param right is the right-hand operand.
     * @return A copy of the current object.
     */
    Matrix& operator=(const Matrix& right);

    /**
     * Equality relational operator: returns <code>true</code> iff both sides
     * are equal. Two matrices are equal if they are of the same type, have equal
     * dimensions and equal values.
     * @param right
     * @return <code>true</code> if the two objects are equal.
     */
    bool operator==(const Matrix& right) const;


    /**
     * Prints out a Matrix object to an output stream.
     * @param os Output stream
     * @param obj Matrix object
     */
    friend std::ostream& operator<<(std::ostream& os, const Matrix& obj);



private:
    /* MatrixFactory is allowed to access these private fields! */
    friend class MatrixFactory;

    size_t m_nrows;         /*< Number of rows */
    size_t m_ncols;         /*< Number of columns */
    bool m_transpose;       /*< Whether this matrix is transposed */
    MatrixType m_type;      /*< Matrix type */

    /* For dense matrices: */

    size_t m_dataLength;        /*< Length of data */
    double *m_data = NULL;      /*< Data */

    /* CSparse members */
    cholmod_triplet *m_triplet = NULL;          /*< Sparse triplets */
    cholmod_sparse *m_sparse = NULL;            /*< A sparse matrix */
    cholmod_factor *m_factor = NULL;            /*< Cholesky factor */
    cholmod_dense *m_dense = NULL;              /*< A dense CHOLMOD matrix */

    
    /* SINGLETON CHOLMOD HANDLE */    
    static cholmod_common *ms_singleton;    /**< Singleton instance of cholmod_common */

    /**
     * Instantiates <code>m_sparse</code> from <code>m_triplet</code>
     * using CHOLMOD's <code>cholmod_triplet_to_sparse</code>. Can only be
     * applied to sparse matrices.
     */
    inline void createSparse();
    
    inline void createTriplet();

    /**
     * Initialize the current matrix (allocate memory etc) for a given number of 
     * rows and columns and a given matrix type.
     * @param nrows Number of rows
     * @param ncols Number of column
     * @param matrixType Matrix type
     */
    inline void init(size_t nrows, size_t ncols, MatrixType matrixType);

    /**
     * Check whether a given pair of indexes is within the matrix bounds.
     * @param i row index
     * @param j column 
     * @return <code>true</code> if (i,j) is within the bounds
     */
    inline bool indexWithinBounds(size_t i, size_t j);

    /**
     * Multiply with a matrix when the left-hand side matrix is dense
     * @param right any right-hand side matrix
     * @return the result of the multiplication (this)*(right) as a new Matrix.
     */
    inline Matrix multiplyLeftDense(const Matrix& right) const;

    /**
     * Multiply with a matrix when the left-hand side matrix is diagonal.
     * @param right any right-hand side matrix
     * @return the result of the multiplication (this)*(right) as a new Matrix.
     */
    inline Matrix multiplyLeftDiagonal(const Matrix& right) const;

    /**
     * Multiply with a matrix when the left-hand side matrix is symmetric.
     * @param right any right-hand side matrix
     * @return the result of the multiplication (this)*(right) as a new Matrix.
     */
    inline Matrix multiplyLeftSymmetric(const Matrix& right) const;

    /**
     * Multiply with a matrix when the left-hand side matrix is sparse.
     * @param right any right-hand side matrix
     * @return the result of the multiplication (this)*(right) as a new Matrix.
     */
    inline Matrix multiplyLeftSparse(Matrix& right);


    /**
     * Custom implementation of matrix-matrix multiplication.
     * @param right RHS
     * @param result the result
     */
    void domm(const Matrix &right, Matrix &result) const;

    /**
     * Storage types for sparse matrix data.
     * This is a private enumeration.
     */
    enum SparseMatrixStorageType {
        CHOLMOD_TYPE_TRIPLET = 444,
        CHOLMOD_TYPE_SPARSE = 555,
        CHOLMOD_TYPE_DENSE = 666,
        CHOLMOD_TYPE_FACTOR = 777
    };

    /**
     * When the matrix is <code>MATRIX_SPARSE</code> this field points to 
     * a CHOLMOD implementation (e.g., <code>cholmod_sparse</code> or
     * <code>cholmod_factor</code>).
     */
    SparseMatrixStorageType m_sparseStorageType;
        
    /**
     * 
     * @param x
     * @return 
     */
    inline double quadFromTriplet(const Matrix& x) const;
    
};

#endif	/* MATRIX_H */

