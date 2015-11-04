/*
 * File:   Matrix.h
 * Author: Pantelis Sopasakis
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
#include <assert.h>

#include "cholmod.h"

#ifdef USE_LIBS
#include <cblas.h>
#include <lapacke.h>
#endif

/**
 * \class Matrix
 * \version version 0.3
 * \ingroup Matrix
 * \date Created on July 7, 2015, 8:02 PM
 *
 * \brief A generic matrix API
 *
 * \note Operations with dense matrices are supported by BLAS and LAPACK routines,
 * while sparse matrix operations are delegated to SuiteSparse routines.
 *
 * <p>A generic matrix which can be an unstructured dense, a structured dense (e.g.,
 * symmetric or lower triangular, stored in packed form) or a dense matrix.
 * This class provides a uniform access framework (an API) to matrix-matrix
 * operations (e.g., addition and multiplication), factorizations and other useful
 * operations.</p>
 *
 * <p>To construct a Matrix you can use one of this class's constructors. However,
 * for sparse matrices it is advisable to use the factory class <code>MatrixFactory</code>.</p>
 *
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

    /**
     * Static method used to destroy the singleton <code>cholmod_handle</code>,
     * if any.
     *
     *
     * @return status this call returns <code>0</code> when it succeeds. See the
     * CHOLMOD documentation for the interpretation of error codes.
     */
    static int destroy_handle();

    /**
     * Types of matrices.
     */
    enum MatrixType {
        MATRIX_DENSE, /**< A dense matrix */
        MATRIX_SPARSE, /**< A sparse matrix (powered by SuiteSparse) */
        MATRIX_DIAGONAL, /**< A diagonal matrix */
        MATRIX_LOWERTR, /**< A lower-triangular matrix */
        MATRIX_SYMMETRIC /**< A symmetric matrix */
    };

    /**
     * Sparse type of matrix
     */
    enum SparseMatrixType {
        SPARSE_UNSYMMETRIC = 0, /**< Not symmetric sparse */
        /**
         * Symmetric sparse and entries in the lower triangular
         * part are transposed and added to the upper triangular part when
         * the matrix is converted to cholmod_sparse form.
         */
        SPARSE_SYMMETRIC_L = 1,
        /**
         * Symmetric sparse and entries in the upper triangular
         * part are transposed and added to the lower triangular part when
         * the matrix is converted to cholmod_sparse form.
         */
        SPARSE_SYMMETRIC_R = -1
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
     *
     * Sets the row and column dimension of this matrix to <code>-1</code> and
     * deletes the matrix data.
     */
    virtual ~Matrix();



    /**
     * Transpose this matrix.
     *
     * Internally, all that this method does is to change the state of this
     * matrix from transposed to non-transposed and vice-versa. Therefore, transposing
     * a matrix using this method does not incur any computation cost.
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
     * Updates a value of the matrix at position <code>(i,j)</code>.
     * 
     * Note that when setting the value <code>A(i,j)</code> of a symmetric matrix,
     * then both <code>A(i,j)</code> and <code>A(j,i)</code> will be set to the 
     * same value.
     * 
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

    std::string getTypeString() const;

    /* Utilities */

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


    bool isSymmetric();

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
     * Computes the quadratic form <code>0.5 * x'*Q*x</code>.
     *
     * <p>Here x is a given vector as <code>Matrix</code>, where Q is the current
     * instance of <code>Matrix</code>.</p>
     *
     * <p>This method can only be applied on square matrices Q while x and q
     * need to be of compatible dimensions.</p>
     *
     * @param x The vector x.
     * @return Scalar 0.5*x'*Q*x as <code>double</code> as a double.
     *
     * \todo Implement quadFromSparse in case the matrix is sparse, but there
     * are no triplets available.
     */
    double quad(Matrix& x);

    /**
     * Computes the quadratic form <code>0.5*x'*Q*x + q'*x</code>.
     *
     * <p>Here x is a given vector as <code>Matrix</code>, where Q is the current
     * instance of <code>Matrix</code>.</p>
     *
     * <p>This method can only be applied on square matrices Q while x and q
     * need to be of combatible dimensions.</p>
     *
     * @param x The vector x.
     * @param q The parameter vector x.
     * @return The result of operation <code>0.5*x'*Q*x + q'*x</code> as a double.
     */
    double quad(Matrix& x, Matrix& q);

    
    /**
     * Converts a vector into a diagonal matrix. This method only applies to 
     * vectors and will throw a <code>std::invalid_argument</code> exception 
     * otherwise.
     */
    void diag();


    /* Operators */

    /**
     * Direct access to the matrix data.
     * Shorthand for <code>matrix.getData()[]</code>; it is however safer to access
     * the matrix entries using <code>get</code> and <code>set</code>.
     *
     * @param sub index
     * @return reference to matrix data
     */
    double &operator[](const size_t sub) const; //overloading []

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
     * Operator *= with a scalar - scale a matrix.
     * @param right
     * @return
     */
    friend Matrix& operator*=(Matrix& obj, double alpha);

    friend Matrix operator*(double alpha, Matrix& obj);

    /**
     * Assignment operator.
     * @param right is the right-hand operand.
     * @return A copy of the current object.
     */
    Matrix& operator=(const double& v);

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

    /**
     * Copy a sub-matrix of the current Matrix object into a new matrix.
     *
     * For dense matrices, this method makes use of LAPACK's <code>dlacpy</code> and
     * is applicable to transpose matrices as well.
     *
     * This method is not supported yet for non-dense matrices (e.g., sparse,
     * diagonal, etc).
     *
     * @param row_start row start
     * @param row_end row end
     * @param col_start column start
     * @param col_end column end
     *
     * @return New matrix object with copied data
     */
    Matrix submatrixCopy(size_t row_start, size_t row_end, size_t col_start, size_t col_end);

    /**
     *
     * Multiply a submatrix of the current <code>Matrix</code> object with a
     * submatrix of a given Matrix.
     *
     * @param rhs Matrix which contains the right-hand side of the multiplication as
     * a submatrix.
     *
     * @param left_row_start
     * @param left_row_end
     * @param left_col_start
     * @param left_col_end
     * @param right_row_start
     * @param right_row_end
     * @param right_col_start
     * @param right_col_end
     * @return The result as a newly created Matrix object.
     */
    Matrix multiplySubmatrix(
            Matrix& rhs,
            const size_t left_row_start,
            const size_t left_row_end,
            const size_t left_col_start,
            const size_t left_col_end,
            const size_t right_row_start,
            const size_t right_row_end,
            const size_t right_col_start,
            const size_t right_col_end);



private:
    /* MatrixFactory is allowed to access these private fields! */
    friend class MatrixFactory;
    friend class CholeskyFactorization;
    friend class LDLFactorization;
    friend class MatrixWriter;

    size_t m_nrows;     /**< Number of rows */
    size_t m_ncols;     /**< Number of columns */
    bool m_transpose;   /**< Whether this matrix is transposed */
    MatrixType m_type;  /**< Matrix type */

    /* For dense matrices: */

    size_t m_dataLength;    /**< Length of data */
    double *m_data;         /**< Data */

    /* CSparse members */
    cholmod_triplet *m_triplet; /**< Sparse triplets */
    cholmod_sparse *m_sparse;   /**< A sparse matrix */
    cholmod_dense *m_dense;     /**< A dense CHOLMOD matrix */


    /* SINGLETON CHOLMOD HANDLE */
    static cholmod_common *ms_singleton; /**< Singleton instance of cholmod_common */

    /**
     * Instantiates <code>m_sparse</code> from <code>m_triplet</code>
     * using CHOLMOD's <code>cholmod_triplet_to_sparse</code>. Can only be
     * applied to sparse matrices.
     */
    void _createSparse();


    /**
     * Creates m_triplet from other existing sparse matrix representations.
     */
    inline void _createTriplet();

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
        CHOLMOD_TYPE_DENSE = 666
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


    inline void _addIJ(size_t i, size_t j, double a);

    /**
     * Adds two vectors as in <code>v1 += v2</code> where <code>v1</code> is a vector
     * to be updated whereas <code>v2</code> is a constant vector.
     * @param len Length of vectors
     * @param pV1 pointer to <code>v1</code>
     * @param pV2 pointer to <code>v2</code>
     */
    static inline void vectorAdd(size_t len, double * pV1, const double * pV2);


    /**
     * Apply operator += on a dense matrix
     * @param rhs right hand side
     */
    inline void _addD(Matrix& rhs); // DENSE + ?

    /**
     * Apply operator += on a symmetric matrix
     * @param rhs right hand side
     */
    inline void _addH(Matrix& rhs); // SYMMETRIC + ?

    /**
     * Apply operator += on a diagonal matrix
     * @param rhs right hand side
     */
    inline void _addX(Matrix& rhs);

    /**
     * Apply operator += on a lower or upper triangular matrix
     * @param rhs right hand side
     */
    inline void _addL(Matrix& rhs);

    /**
     * Apply operator += on a sparse matrix
     * @param rhs right hand side
     */
    inline void _addS(Matrix& rhs);


};

#endif	/* MATRIX_H */
