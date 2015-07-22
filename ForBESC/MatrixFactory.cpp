/* 
 * File:   MatrixFactory.cpp
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

#include <set>

#include "MatrixFactory.h"
#include "Matrix.h"


Matrix MatrixFactory::MakeIdentity(int n, float alpha) {
    Matrix mat(n, n, Matrix::MATRIX_DIAGONAL);
    for (int i = 0; i < n; i++) {
        mat[i] = alpha;
    }
    return mat;
}

typedef std::pair<int, int> nice_pair;

Matrix MatrixFactory::MakeRandomSparse(int nrows, int ncols, int nnz, float offset, float scale) {
    Matrix R = MakeSparse(nrows, ncols, nnz, Matrix::SPARSE_UNSYMMETRIC);
    std::srand(unsigned ( std::time(0)));
    std::set<nice_pair> s;
    nice_pair p;
    while (true) { // construct pairs
        p.first = (std::rand() % (int) (nrows));
        p.second = (std::rand() % (int) (ncols));
        s.insert(p);
        if (s.size() == nnz) {
            break;
        }
    }
    float rand;
    for (std::set<nice_pair>::iterator it = s.begin(); it != s.end(); it++) {
        rand = offset + scale * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        R.set(it->first, it->second, rand);
    }
    return R;
}

Matrix MatrixFactory::MakeRandomMatrix(int nrows, int ncols, float offset, float scale, Matrix::MatrixType type) {
    int len = 0;
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
            
            break;
    }
    Matrix mat(nrows, ncols, type);
    for (int j = 0; j < len; j++) {
        mat[j] = offset + scale * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
    }
    return mat;
}


Matrix MatrixFactory::MakeSparse(int nrows, int ncols, int max_nnz, Matrix::SparseMatrixType stype) {
    Matrix matrix(nrows, ncols, Matrix::MATRIX_SPARSE);
    matrix.m_triplet = cholmod_allocate_triplet(nrows, ncols, max_nnz, stype, CHOLMOD_REAL, Matrix::cholmod_handle());
    matrix.m_sparseStorageType = Matrix::CHOLMOD_TYPE_TRIPLET;
    return matrix;
}

Matrix MatrixFactory::MakeSparseSymmetric(int nrows, int ncols, int max_nnz) {
    return MakeSparse(nrows, ncols, max_nnz, Matrix::SPARSE_SYMMETRIC_L);
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






