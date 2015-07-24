/* 
 * File:   MatrixFactory.h
 * Author: chung
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
#include <vector>       // std::vector
#include <algorithm>    // std::random_shuffle
#include <ctime>        // std::time
#include <functional>

class MatrixFactory {
public:
    
    cholmod_common COMMON;
    

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
     * Constructs an identity matrix of type <code>Matrix::MATRIX_DIAGONAL</code>.
     * 
     * @param n number of diagonal elements
     * @param alpha scaling factor
     * @return Diagonal matrix of the form <code>alpha*I</code>, where <code>alpha</code>
     * is a given scalar.
     */
    static Matrix MakeIdentity(size_t n, float alpha);
        
    
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

private:

};

#endif	/* MATRIXFACTORY_H */

