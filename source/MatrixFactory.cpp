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

Matrix MatrixFactory::MakeIdentity(size_t n, float alpha) {
    Matrix mat(n, n, Matrix::MATRIX_DIAGONAL);
    for (size_t i = 0; i < n; i++) {
        mat[i] = alpha;
    }
    return mat;
}

Matrix MatrixFactory::MakeRandomSparse(size_t nrows, size_t ncols, size_t nnz, float offset, float scale) {
    std::srand((unsigned int) std::time(0));
    Matrix R = MakeSparse(nrows, ncols, nnz, Matrix::SPARSE_UNSYMMETRIC);
    std::set<nice_pair> s;
    nice_pair p;
    while (true) { // construct pairs
        p.first = (std::rand() %  nrows);
        p.second = (std::rand() %  ncols);
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
    std::srand((unsigned int) std::time(0));
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
    }
    Matrix mat(nrows, ncols, type);
    for (size_t j = 0; j < len; j++) {
        mat[j] = offset + scale * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
    }
    return mat;
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