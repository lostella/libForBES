/* 
 * File:   MatrixFactory.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 12, 2015, 7:50 PM
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

#ifndef MATRIXFACTORY_H
#define	MATRIXFACTORY_H

#include "Matrix.h"

/**
 * \class MatrixFactory
 * \brief Use this class to create instances of %Matrix
 * \version version 0.1
 * \ingroup Matrix-group
 * \date Created on July 12, 2015, 7:50 PM
 * \author <a href="http://dysco.imtlucca.it/sopasakis">Pantelis Sopasakis</a>
 * \example shallow_example.cpp
 * 
 * 
 * \todo Create matrix constructor from double* using a pointer so as not to allocate
 * any memory inside the Matrix object.
 * 
 */
class MatrixFactory {
public:


    /**
     * Returns a random matrix as an instance of <code>Matrix</code>.
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @param offset Random data follow a statistical distribution of the form 
     * <code>offset + scaling_factor * r</code>, where <code>r</code> is a random
     * variable in [0,1].     
     * @param scale The aforementioned scaling factor.
     * @param type matrix type.
     * @return A random matrix object.
     */
    static Matrix MakeRandomMatrix(size_t nrows, size_t ncols, float offset, float scale, Matrix::MatrixType type);


    /**
     * Returns a random matrix as an instance of <code>Matrix</code> with type
     * \link Matrix::MATRIX_DENSE MATRIX_DENSE\endlink..
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @param offset Random data follow a statistical distribution of the form 
     * <code>offset + scaling_factor * r</code>, where <code>r</code> is a random
     * variable in [0,1].     
     * @param scale The aforementioned scaling factor.
     * @return A random dense matrix object.
     */
    static Matrix MakeRandomMatrix(size_t nrows, size_t ncols, float offset, float scale);


    /**
     * Constructs an identity matrix of type <code>Matrix::MATRIX_DIAGONAL</code>.
     * 
     * @param n number of diagonal elements
     * @param alpha scaling factor
     * @return Diagonal matrix of the form <code>alpha*I</code>, where <code>alpha</code>
     * is a given scalar.
     */
    static Matrix MakeIdentity(size_t n, double alpha);


    /**
     * Creates a sparse matrix of given dimensions, maximum number of non-zero 
     * elements and sparsity pattern (symmetric or not).
     * @param nrows number of rows
     * @param ncols number of columns.
     * @param max_nnz maximum number of non-zero elements
     * @param stype symmetry type
     * @return Allocated sparse matrix
     */
    static Matrix MakeSparse(size_t nrows, size_t ncols, size_t max_nnz, Matrix::SparseMatrixType stype);

    /**
     * Creates a sparse symmetric matrix of given dimensions.
     * @param n matrix size (matrix is square)
     * @param max_nnz maximum number of non-zero elements
     * @return Allocated sparse matrix
     */
    static Matrix MakeSparseSymmetric(size_t n, size_t max_nnz);

    /**
     * Allocates a sparse matrix of given dimensions and instantiates it with
     * random entries at random positions. The client needs to specify the
     * number of non-zero elements of the matrix.
     * 
     * Random data follow a statistical distribution of the form 
     * <code>offset + scale * r</code>, where <code>r</code> is a random
     * variable in [0,1]. Parameters <code>offset</code> and <code>scale</code> 
     * are provided as input arguments.
     * 
     * The allocated matrix assumes that its <em>lower-triangular part</em> is 
     * provided only. 
     * 
     * @param nrows number of rows
     * @param ncols number of columns.
     * @param nnz number of non-zero elements (must be less than <code>nrows*ncols</code>)
     * @param offset random number offset
     * @param scale random number scaling factor
     * @return Allocated random sparse matrix
     */
    static Matrix MakeRandomSparse(size_t nrows, size_t ncols, size_t nnz, float offset, float scale);

    /**
     * Reads a sparse matrix from a file and returns a <code>Matrix</code> object.
     * 
     * The caller provides a pointer to a <code>FILE</code> which can be 
     * created using <code>fopen</code>. The matrix file is an ASCII text file
     * with space-separated values (triplets). The first triplet defines the
     * matrix dimensions as <code>(n_rows, n_columns, n_non_zero)</code>.
     * Subsequent triplets are in the form <code>(i, j, value)</code> where
     * i: [0,..., n_rows - 1] and j: [0, ..., n_columns - 1].
     * 
     * @param fp A pointer to a file object
     * @return Sparse matrix constructed from a file
     * @throw cholmod_error an exception is thrown when the file format is wrong.
     */
    static Matrix ReadSparse(FILE *fp);

    /**
     * Creates a shallow matrix given a %Matrix object.
     * 
     * @param mat original matrix object
     * @return  shallow matrix
     */
    static Matrix ShallowMatrix(const Matrix& mat);
    
    /**
     * Creates a <em>shallow vector</em> from a given %Matrix object.
     * Shallow vectors do not allocate space for their data, but instead point to 
     * an allocated space in memory.
     * 
     * \note Since the internal data of a shallow vector point to another object's
     * data - here to the data of <code>vector</code>. Any modification in the original
     * (actual) data will be visible from the shallow copy. However, reshape operations
     * and transposition of the original matrix will not affect its shallow copy.
     * 
     * \exception std::invalid_argument in case the provided %Matrix object <code>vector</code>
     * is not a dense matrix. Note that it is allowed for the provided %Matrix object to 
     * be either a column or a row vector.
     * 
     * \post
     * The %Matrix object that is returned by this method "points" (internally) to the
     * original %Matrix object; therefore, any modification of the new object's internal state
     * affects directly the state of the original object.
     * 
     * @param vector any dense vector (either column or row vector)
     * @param size size/length of the shallow vector to be created
     * @param offset offset with respect to the original vector
     * @return shallow vector
     */
    static Matrix ShallowVector(const Matrix& vector, size_t size, size_t offset);

    /**
     * Creates a <em>shallow vector</em> from a given %Matrix object.
     * Shallow vectors do not allocate space for their data, but instead point to 
     * an allocated space in memory.
     * 
     * \note Since the internal data of a shallow vector point to another object's
     * data - here to the data of <code>vector</code>. Any modification in the original
     * (actual) data will be visible from the shallow copy. However, reshape operations
     * and transposition of the original matrix will not affect its shallow copy.
     * 
     * \exception std::invalid_argument in case the provided %Matrix object <code>vector</code>
     * is not a dense matrix. Note that it is allowed for the provided %Matrix object to 
     * be either a column or a row vector.
     * 
     * \post
     * The %Matrix object that is returned by this method "points" (internally) to the
     * original %Matrix object; therefore, any modification of the new object's internal state
     * affects directly the state of the original object.
     * 
     * @param vector any dense vector (either column or row vector)
     * @param offset offset with respect to the original vector
     * @return shallow vector
     * 
     * \sa \link Matrix::submatrixCopy submatrixCopy\endlink
     */
    static Matrix ShallowVector(const Matrix& vector, size_t offset);

    /**
     * Creates a <em>shallow vector</em> from a given pointer-to-double.
     * Shallow vectors do not allocate space for their data, but instead point to 
     * an allocated space in memory.
     * 
     * \post
     * The %Matrix object that is returned by this method "points" (internally) to the
     * original %Matrix object; therefore, any modification of the new object's internal state
     * affects directly the state of the original object.
     * 
     * @param data pointer to data
     * 
     * @param size size of the created shallow vector. This should be of course
     * smaller than or equal to the size of the original data.
     * 
     * @param offset offset with respect to the original pointer
     * 
     * @return shallow vector
     * 
     * \sa \link Matrix::submatrixCopy submatrixCopy\endlink
     */
    static Matrix ShallowVector(double * data, size_t size, size_t offset);
    
    /**
     * Creates an empty shallow vector/matrix which does not point anywhere. 
     * The size of this %Matrix object will be 0-by-0.
     * 
     * @return empty shallow vector
     */
    static Matrix ShallowVector();


private:

};

#endif	/* MATRIXFACTORY_H */

