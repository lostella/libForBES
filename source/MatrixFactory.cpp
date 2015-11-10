/* 
 * File:   MatrixFactory.cpp
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

#include <set>

#include "MatrixFactory.h"
#include "Matrix.h"

typedef std::pair<size_t, size_t> nice_pair;

Matrix MatrixFactory::MakeIdentity(size_t n, double alpha) {
    Matrix mat(n, n, Matrix::MATRIX_DIAGONAL);
    for (size_t i = 0; i < n; i++) {
        mat[i] = alpha;
    }
    return mat;
}

Matrix MatrixFactory::MakeRandomSparse(size_t nrows, size_t ncols, size_t nnz, float offset, float scale) {
    std::srand(std::time(0));
    Matrix R = MakeSparse(nrows, ncols, nnz, Matrix::SPARSE_UNSYMMETRIC);
    std::set<nice_pair> s;
    nice_pair p;
    while (true) { // construct pairs
        p.first = (std::rand() % nrows);
        p.second = (std::rand() % ncols);
        s.insert(p);
        if (s.size() == nnz) {
            break;
        }
    }
    for (std::set<nice_pair>::iterator it = s.begin(); it != s.end(); ++it) {
        float rand;
        rand = offset + scale * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        R.set(it->first, it->second, rand);
    }
    return R;
}

Matrix MatrixFactory::MakeRandomMatrix(size_t nrows, size_t ncols, float offset, float scale, Matrix::MatrixType type) {
    std::srand(std::time(0));
    size_t len = 0;
    switch (type) {
        case Matrix::MATRIX_DENSE:
            len = nrows * ncols;
            break;
        case Matrix::MATRIX_LOWERTR:
        case Matrix::MATRIX_SYMMETRIC:
            len = nrows * (nrows + 1) / 2;
            break;
        case Matrix::MATRIX_DIAGONAL:
            len = nrows;
            break;
        case Matrix::MATRIX_SPARSE:
            return MakeRandomSparse(nrows, ncols, std::ceil((nrows * ncols) / 3), offset, scale);
            break;
        default:
            throw std::logic_error("unsupported");
    }
    Matrix mat(nrows, ncols, type);
    for (size_t j = 0; j < len; j++) {
        mat[j] = offset + scale * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
    }
    return mat;
}

Matrix MatrixFactory::MakeRandomMatrix(size_t nrows, size_t ncols, float offset, float scale) {
    return MatrixFactory::MakeRandomMatrix(nrows, ncols, offset, scale, Matrix::MATRIX_DENSE);
}

Matrix MatrixFactory::MakeSparse(size_t nrows, size_t ncols, size_t max_nnz, Matrix::SparseMatrixType stype) {
    Matrix matrix(nrows, ncols, Matrix::MATRIX_SPARSE);
    matrix.m_triplet = cholmod_allocate_triplet(nrows, ncols, max_nnz, stype, CHOLMOD_REAL, Matrix::cholmod_handle());
    matrix.m_sparseStorageType = Matrix::CHOLMOD_TYPE_TRIPLET;
    return matrix;
}

Matrix MatrixFactory::MakeSparseSymmetric(size_t n, size_t max_nnz) {
    return MakeSparse(n, n, max_nnz, Matrix::SPARSE_SYMMETRIC_L);
}

Matrix MatrixFactory::ReadSparse(FILE* fp) {
    cholmod_sparse *sp;
    sp = cholmod_read_sparse(fp, Matrix::cholmod_handle());

    Matrix mat(sp->nrow, sp->ncol, Matrix::MATRIX_SPARSE);
    mat.m_sparse = sp;
    mat.m_sparseStorageType = Matrix::CHOLMOD_TYPE_SPARSE;
    mat.m_triplet = cholmod_sparse_to_triplet(sp, Matrix::cholmod_handle());
    return mat;
}


Matrix MatrixFactory::ShallowMatrix(const Matrix& orig) {   
    Matrix mat_shallow = Matrix(true);
    mat_shallow.m_transpose = orig.m_transpose;
    mat_shallow.m_nrows = orig.m_nrows;
    mat_shallow.m_ncols = orig.m_ncols;
    mat_shallow.m_dataLength = orig.m_dataLength;
    mat_shallow.m_delete_data = false;
    mat_shallow.m_data = orig.m_data;
    mat_shallow.m_type = orig.m_type;
    mat_shallow.m_sparseStorageType = orig.m_sparseStorageType;
    mat_shallow.m_dense = orig.m_dense;
    mat_shallow.m_triplet = orig.m_triplet;
    mat_shallow.m_sparse = orig.m_sparse;
    return mat_shallow;
}

Matrix MatrixFactory::ShallowVector(const Matrix& orig, size_t size, size_t offset) {
    if (orig.getType() != Matrix::MATRIX_DENSE) {
        throw std::invalid_argument("This method can only be applied to dense vectors");
    }
    if (orig.getNcols() != 1 && orig.getNrows() != 1) {
        throw std::invalid_argument("This method can only be applied to vectors");
    }
    Matrix v_shallow = Matrix(true);
    v_shallow.m_transpose = orig.m_transpose;
    v_shallow.m_nrows = orig.m_transpose ? 1 : size;
    v_shallow.m_ncols = orig.m_transpose ? size : 1;
    v_shallow.m_dataLength = v_shallow.m_nrows * v_shallow.m_ncols;
    v_shallow.m_delete_data = false;
    v_shallow.m_data = orig.m_data + offset;
    return v_shallow;
}

Matrix MatrixFactory::ShallowVector(const Matrix& orig, size_t offset) {
    return MatrixFactory::ShallowVector(orig, orig.length() - offset, offset);
}

Matrix MatrixFactory::ShallowVector(double* data, size_t size, size_t offset) {
    Matrix v_shallow = Matrix(true);
    v_shallow.m_transpose = false;
    v_shallow.m_nrows = size;
    v_shallow.m_ncols = 1;
    v_shallow.m_dataLength = v_shallow.m_nrows;
    v_shallow.m_delete_data = false;
    v_shallow.m_data = data + offset;
    return v_shallow;
}

Matrix MatrixFactory::ShallowVector() {
    Matrix v_shallow = Matrix(true);
    return v_shallow;
}




