/* 
 * File:   Matrix.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on June 30, 2015, 12:34 PM
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

#include "Matrix.h"
#include <iostream>
#include <stdexcept>
#include <complex>
#include <cstdlib>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <assert.h>
#include <limits>

#ifdef USE_LIBS
#include <cblas.h>
#include <lapacke.h>
#endif

/* STATIC MEMBERS */

cholmod_common* Matrix::ms_singleton = NULL;

cholmod_common* Matrix::cholmod_handle() {
    if (ms_singleton == NULL) {
        ms_singleton = new cholmod_common;
        cholmod_start(ms_singleton);
    }
    return ms_singleton;
}

int Matrix::destroy_handle() {
    if (ms_singleton == NULL) {
        return 0;
    }
    int status = cholmod_finish(ms_singleton);
    ms_singleton = NULL;
    return status;
}

/********* CONSTRUCTORS ************/
Matrix::Matrix() {
    m_nrows = 0;
    m_ncols = 0;
    m_data = new double[1];
    *m_data = 0;
    m_type = MATRIX_DENSE;
    m_dataLength = 0;
    m_transpose = false;
    m_triplet = NULL;
    m_sparse = NULL;
    m_dense = NULL;
    m_sparseStorageType = CHOLMOD_TYPE_TRIPLET;
    m_delete_data = true;
}

Matrix::Matrix(std::pair<size_t,size_t> dimensions){
    init(dimensions.first, dimensions.second, MATRIX_DENSE);
}

Matrix::Matrix(size_t nr, size_t nc) {
    init(nr, nc, MATRIX_DENSE);
}

Matrix::Matrix(size_t nr, size_t nc, MatrixType mType) {
    init(nr, nc, mType);
}

Matrix::Matrix(size_t nr, size_t nc, const double * dat) {
    init(nr, nc, MATRIX_DENSE);
    for (size_t j = 0; j < nc * nr; j++) {
        m_data[j] = dat[j];
    }
}

Matrix::Matrix(size_t nr, size_t nc, const double * dat, MatrixType mType) {
    init(nr, nc, mType);
    for (size_t j = 0; j < length(); j++) {
        m_data[j] = dat[j];
    }
}

Matrix::Matrix(const Matrix& orig) {
    m_ncols = orig.m_ncols;
    m_nrows = orig.m_nrows;
    m_transpose = orig.m_transpose;
    m_data = NULL;
    m_delete_data = orig.m_delete_data;
    m_triplet = NULL;
    m_sparse = NULL;
    m_dense = NULL;
    m_type = orig.m_type;
    if (orig.m_type != MATRIX_SPARSE) {
        size_t n = orig.m_dataLength;
        if (n == 0) {
            n = 1;
        }
        m_data = new double[n];        
        memcpy(m_data, orig.m_data, n*sizeof(double));
        m_dataLength = orig.m_dataLength;
        m_delete_data = true;
    } else {
        if (orig.m_triplet != NULL) {
            m_triplet = cholmod_copy_triplet(orig.m_triplet, Matrix::cholmod_handle());
        }
        if (orig.m_sparse != NULL) {
            m_sparse = cholmod_copy_sparse(orig.m_sparse, Matrix::cholmod_handle());
        }
        if (orig.m_dense != NULL) {
            m_dense = cholmod_copy_dense(orig.m_dense, Matrix::cholmod_handle());
        }
    }
    m_sparseStorageType = orig.m_sparseStorageType;
}

/********* DENSTRUCTOR ************/
Matrix::~Matrix() {
    m_ncols = 0;
    m_nrows = 0;
    if (m_data != NULL && m_delete_data) {
        delete[] m_data;
    }
    m_data = NULL; /* so that we don't double-free */
    m_delete_data = false; /* for extra safety (just in case) */
    if (m_triplet != NULL) {
        cholmod_free_triplet(&m_triplet, Matrix::cholmod_handle());
        m_triplet = NULL;
    }
    if (m_sparse != NULL) {
        cholmod_free_sparse(&m_sparse, Matrix::cholmod_handle());
        m_sparse = NULL;
    }
    if (m_dense != NULL) {
        m_dense->x = NULL;
        cholmod_free_dense(&m_dense, Matrix::cholmod_handle());
        m_dense = NULL;
    }
}

/********* GETTERS/SETTERS ************/
size_t Matrix::getNcols() const {
    return m_ncols;
}

size_t Matrix::getNrows() const {
    return m_nrows;
}

double * Matrix::getData() {
    return m_data;
}

bool Matrix::isEmpty() const {
    return (m_nrows == 0 || m_ncols == 0);
}

bool Matrix::isColumnVector() const {
    return this -> m_ncols == 1;
}

bool Matrix::isRowVector() const {
    return this -> m_nrows == 1;
}

size_t Matrix::length() const {
    return m_dataLength;
}

/********* OTHER METHODS ************/

void Matrix::transpose() {
    if (m_type == MATRIX_DIAGONAL || m_type == MATRIX_SYMMETRIC) {
        return;
    }

    if (m_type == MATRIX_SPARSE) {
        _createSparse();
        m_sparse = cholmod_transpose(m_sparse, 1, cholmod_handle());
    }
    if (this -> m_transpose) {
        this -> m_transpose = false;
    } else {
        this -> m_transpose = true;
    }
    std::swap(this -> m_ncols, this -> m_nrows);
}

int Matrix::reshape(size_t nrows, size_t ncols) {
    size_t new_size = nrows * ncols;
    if (new_size == 0) {
        return -1;
    }
    if (new_size > length()) {
        return -2;
    }
    this -> m_nrows = nrows;
    this -> m_ncols = ncols;
    return 0;
}

double Matrix::get(const size_t i) const {
    return m_data[i];
}

double Matrix::get(const size_t i, const size_t j) const {
    //LCOV_EXCL_START
    if (isEmpty()) {
        throw std::out_of_range("Method get(size_t, size_t) applied to an empty matrix");
    }
    if (i >= getNrows() || j >= getNcols()) {
        throw std::out_of_range("Index out of range!");
    }
    //LCOV_EXCL_STOP
    if (m_type == MATRIX_DENSE) {
        return !m_transpose ? m_data[i + j * m_nrows] : m_data[j + i * m_ncols];
    } else if (m_type == MATRIX_DIAGONAL) {
        if (i == j) {
            return m_data[i];
        } else {
            return 0.0;
        }
    } else if (m_type == MATRIX_SYMMETRIC) {
        size_t i_ = std::max(i, j);
        size_t j_ = std::min(i, j);
        return m_data[i_ + m_nrows * j_ - j_ * (j_ + 1) / 2];
    } else if (m_type == MATRIX_LOWERTR) {
        return m_transpose
                ? (j >= i) ? m_data[j + m_ncols * i - i * (i + 1) / 2] : 0.0
                : (i >= j) ? m_data[i + m_nrows * j - j * (j + 1) / 2] : 0.0;
    } else {
        /* if (m_type == MATRIX_SPARSE) */
        //LCOV_EXCL_START
        if (m_triplet == NULL) {
            throw std::logic_error("not supported yet");
        }
        //LCOV_EXCL_STOP
        double val = 0.0;
        int i_ = m_transpose ? j : i;
        int j_ = m_transpose ? i : j;
        for (size_t k = 0; k < m_triplet->nnz; k++) {
            if (i_ == (static_cast<int*> (m_triplet->i))[k]) {
                if (j_ == (static_cast<int*> (m_triplet->j))[k]) {
                    val = (static_cast<double*> (m_triplet->x))[k];
                    break;
                }
            }
        }
        return val;
    }

} /* END GET */

void Matrix::set(size_t i, size_t j, double v) {
    //LCOV_EXCL_START
    if (!indexWithinBounds(i, j)) {
        throw std::out_of_range("Index out of range!");
    }
    //LCOV_EXCL_STOP
    if (m_type == MATRIX_DENSE) {
        if (m_transpose) {
            m_data[j + i * m_ncols] = v;
        } else {
            m_data[i + j * m_nrows] = v;
        }
    } else if (m_type == MATRIX_DIAGONAL && i == j) {
        m_data[i] = v;
    } else if (m_type == MATRIX_LOWERTR) {
        m_data[i + m_nrows * j - j * (j + 1) / 2] = v;
    } else if (m_type == MATRIX_SYMMETRIC) { /* sets A(i,j) = A(j,i) = v */
        int i_ = std::max(i, j);
        int j_ = std::min(i, j);
        m_data[i_ + m_nrows * j_ - j_ * (j_ + 1) / 2] = v;
    } else if (m_type == MATRIX_SPARSE) {
        if (m_triplet == NULL) { /* Create triplets if they don't exist */
            _createTriplet();
        }
        if (m_triplet->nnz == m_triplet->nzmax) { /* max NNZ exceeded */
            cholmod_reallocate_triplet(m_triplet->nzmax + 1, m_triplet, Matrix::cholmod_handle());
        }

        int k_found = -1;
        for (size_t s = 0; s < m_triplet->nnz; s++) {
            if (i == (static_cast<int*> (m_triplet->i))[s] && j == (static_cast<int*> (m_triplet->j))[s]) {
                k_found = s;
                break;
            }
        }

        if (k_found == -1) {
            (static_cast<int*> (m_triplet->i))[m_triplet->nnz] = i;
            (static_cast<int*> (m_triplet->j))[m_triplet->nnz] = j;
            (static_cast<double*> (m_triplet->x))[m_triplet->nnz] = v;
            (m_triplet->nnz)++;
        } else {
            (static_cast<int*> (m_triplet->i))[k_found] = i;
            (static_cast<int*> (m_triplet->j))[k_found] = j;
            (static_cast<double*> (m_triplet->x))[k_found] = v;
        }
        /* Invalidate alternative sparse representations */
        if (m_sparse != NULL) {
            cholmod_free_sparse(&m_sparse, Matrix::cholmod_handle());
        }
        if (m_dense != NULL) {
            cholmod_free_dense(&m_dense, Matrix::cholmod_handle());
        }
        m_sparse = NULL;
        m_dense = NULL;
    } else {
        //LCOV_EXCL_START
        throw std::invalid_argument("Illegal operation");
        //LCOV_EXCL_STOP
    }

}

double Matrix::quadFromTriplet(const Matrix& x) const {
    double r = 0.0;
    for (size_t k = 0; k < m_triplet->nnz; k++) {
        r += x.get(static_cast<int*> (m_triplet->i)[k], 0) *
                x.get(static_cast<int*> (m_triplet->j)[k], 0) *
                (static_cast<double*> (m_triplet->x))[k];
    }
    return r;
}

double Matrix::quad(Matrix & x) const {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("Method `quadratic` can only be applied to vectors!");
    }
    if (getNcols() != getNrows()) {
        throw std::invalid_argument("Method `quadratic` can only be applied on square matrices!");
    }
    if (x.getNrows() != m_ncols) {
        throw std::invalid_argument("The argument of quad(Matrix&) is not of appropriate dimension.");
    }
    //LCOV_EXCL_STOP
    double result = 0.0;

    if (MATRIX_DENSE == m_type || MATRIX_LOWERTR == m_type) { /* DENSE or LOWER TRIANGULAR */
        for (size_t j = 0; j < m_ncols; j++) {
            for (size_t i = 0; i < m_nrows; i++) {
                result += x[i] * (!m_transpose ? m_data[i + j * m_nrows] : m_data[j + i * m_ncols]) * x[j];
            }
        }
    } else if (MATRIX_DIAGONAL == m_type) { /* DIAGONAL */
        for (size_t i = 0; i < m_nrows; i++) {
            result += x[i] * x[i] * m_data[i];
        }
    } else if (MATRIX_SYMMETRIC == m_type) { /* SYMMETRIC */
        for (size_t i = 0; i < m_nrows; i++) {
            for (size_t j = 0; j < i; j++) {
                result += 2 * x[i] * get(i, j) * x[j];
            }
            result += x[i] * get(i, i) * x[i];
        }
    } else if (MATRIX_SPARSE == m_type) { /* SPARSE */
        if (m_triplet != NULL) {
            result = quadFromTriplet(x);
        } else {
            /* \todo Implement quadFromSparse */
            //LCOV_EXCL_START
            throw std::logic_error("Quad on sparse matrix - no triplets found (not implemented yet)");
            //LCOV_EXCL_STOP
        }
    }
    result /= 2.0;
    return result;
}

double Matrix::quad(Matrix& x, Matrix & q) const {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("Method `quadratic` can only be applied to vectors!");
    }
    if (!q.isColumnVector()) {
        throw std::invalid_argument("Parameter q needs to be a column vector!");
    }
    if (getNcols() != getNrows()) {
        throw std::invalid_argument("Method `quadratic` can only be applied on square matrices!");
    }
    if (x.getNrows() != m_ncols) {
        throw std::invalid_argument("The argument of quad(Matrix&) is not of appropriate dimension.");
    }
    //LCOV_EXCL_STOP
    double t;
    Matrix r;
    r = q*x;
    assert(!r.isEmpty());
    t = r.get(0, 0) + quad(x); /* q*x + (1/2) x'*Q*x */

    return t;
}

void Matrix::plusop() {
    if (m_type != Matrix::MATRIX_SPARSE) {
        for (size_t i = 0; i < length(); i++) {
            if (m_data[i] < 0) {
                m_data[i] = 0.0;
            }
        }
    } else {
        if (m_triplet != NULL) {
            for (size_t k = 0; k < m_triplet->nnz; k++) {
                double * val = (static_cast<double*> (m_triplet->x)) + k;
                if (*val < 0) {
                    *val = 0.0;
                }
            }
        } else if (m_sparse != NULL) {
            for (size_t j = 0; j < m_ncols; j++) {
                int p = (static_cast<int*> (m_sparse->p))[j];
                int pend = (m_sparse->packed == 1)
                        ? ((static_cast<int*> (m_sparse->p))[j + 1])
                        : p + (static_cast<int*> (m_sparse->nz))[j];
                for (; p < pend; p++) {
                    double * val = (static_cast<double*> (m_sparse->x)) + p;
                    if (*val < 0) {
                        *val = 0.0;
                    }
                }
            }
        }
    }
}

void Matrix::plusop(Matrix* mat) {
    if (m_type != Matrix::MATRIX_SPARSE) {
        if (length() != mat->length()) {
            throw std::invalid_argument("Input matrix allocation/size error");
        }
        for (size_t i = 0; i < length(); i++) {
            mat->m_data[i] = m_data[i] < 0.0 ? 0.0 : m_data[i];
        }
    } else {
        throw std::logic_error("Not implemented yet");
    }
}

/********* OPERATORS ************/
bool Matrix::operator==(const Matrix & right) const {
    const double tol = 1e-9;
    bool result = (m_type == right.m_type) &&
            (m_ncols == right.m_ncols) &&
            (m_nrows == right.m_nrows);
    for (unsigned int i = 0; i < m_nrows; i++) {
        for (size_t j = 0; j < m_ncols; j++) {
            result = result && (std::abs(get(i, j) - right.get(i, j)) < tol);
        }
    }
    return result;
}

//LCOV_EXCL_START

std::ostream& operator<<(std::ostream& os, const Matrix & obj) {
    os << "\nMatrix " << obj.m_nrows << "x" << obj.m_ncols << std::endl;
    if (obj.m_transpose) {
        os << "Stored as transpose : YES\n";
    }
    const char * const types[] = {"Dense", "Sparse", "Diagonal", "Lower Triangular", "Symmetric"};
    os << "Type: " << types[obj.m_type] << std::endl;
    if (obj.m_type == Matrix::MATRIX_SPARSE && obj.m_triplet == NULL && obj.m_sparse != NULL) {
        os << "Storage type: Packed Sparse" << std::endl;
        for (size_t j = 0; j < obj.m_ncols; j++) {
            int p = (static_cast<int*> (obj.m_sparse->p))[j];
            int pend = (obj.m_sparse->packed == 1)
                    ? ((static_cast<int*> (obj.m_sparse->p))[j + 1])
                    : p + (static_cast<int*> (obj.m_sparse->nz))[j];
            for (; p < pend; p++) {
                os << "(" << (static_cast<int*> (obj.m_sparse->i))[p] << "," << j << ")  : "
                        << std::setw(9) << std::setprecision(4) << (static_cast<double*> (obj.m_sparse->x))[p] << std::endl;
            }
        }
        return os;
    }
    for (size_t i = 0; i < obj.m_nrows; i++) {
        for (size_t j = 0; j < obj.m_ncols; j++) {
            if (obj.m_type == Matrix::MATRIX_SPARSE && obj.get(i, j) == 0) {
                os << std::setw(8) << "0";
            } else {

                os << std::setw(8) << std::setprecision(4) << obj.get(i, j);
            }
        }
        os << std::endl;
    }
    return os;
}
//LCOV_EXCL_STOP

double &Matrix::operator[](size_t sub) const {
    //LCOV_EXCL_START
//    if (sub >= length()) {
//        throw std::out_of_range("Exception: Index out of range for Matrix");
//    }
    //LCOV_EXCL_STOP
    return m_data[sub];
}

inline void Matrix::_addIJ(size_t i, size_t j, double a) {
    set(i, j, get(i, j) + a); // A(i,j) += a
}

inline void Matrix::_addIJ(size_t i, size_t j, double a, double gamma) {
    set(i, j, gamma * get(i, j) + a); // A(i,j) += a
}

inline void Matrix::vectorAdd(size_t len, double* pV1, const double* pV2) {
#ifdef USE_LIBS 
    cblas_daxpy(len, 1.0, pV2, 1, pV1, 1); // data = data + right.data
#else
    for (int i = 0; i < len; i++)
        pV1[i] += pV2[i];
#endif
}

Matrix& Matrix::operator+=(Matrix & right) {
    //LCOV_EXCL_START
    if (m_ncols != right.m_ncols || m_nrows != right.m_nrows) {
        throw std::invalid_argument("Incompatible dimensions while using +=!");
    }
    //LCOV_EXCL_STOP

    if (&right == this) {
        *this *= 2.0;
        return *this;
    }

    const double alpha = 1.0;
    const double gamma = 1.0;
    add(*this, alpha, right, gamma);
    return *this;
}

Matrix & Matrix::operator-=(Matrix & right) {
    //LCOV_EXCL_START
    if (m_ncols != right.m_ncols || m_nrows != right.m_nrows) {
        throw std::invalid_argument("Incompatible dimensions while using +=!");
    }
    //LCOV_EXCL_STOP

    const double alpha = -1.0;
    const double gamma = 1.0;
    add(*this, alpha, right, gamma);
    return *this;
}

Matrix Matrix::operator+(Matrix & right) const {
    if (this->getNrows() != right.getNrows() || this->getNcols() != right.getNcols()) {
        throw std::invalid_argument("Addition of matrices of incompatible dimensions!");
    }
    Matrix result(*this); // Make a copy of myself.
    result += right;

    return result;
}

Matrix Matrix::operator-(Matrix & right) const {
    if (this->getNrows() != right.getNrows() || this->getNcols() != right.getNcols()) {
        throw std::invalid_argument("Addition of matrices of incompatible dimensions!");
    }
    Matrix result(*this); // Make a copy of myself.      
    result -= right;

    return result;
}

Matrix Matrix::operator*(Matrix & right) {
    if (!(getType() == Matrix::MATRIX_SPARSE && right.getType() == Matrix::MATRIX_SPARSE) &&
            isColumnVector() && right.isColumnVector() && length() == right.length()) {
        double t = 0.0;
        // multiplication of two column vectors = dot product
        Matrix r(1, 1);
        for (size_t i = 0; i < m_nrows * m_ncols; i++) {
            t += m_data[i] * right.m_data[i];
        }
        r[0] = t;
        return r;
    }
    if (!(isColumnVector() && right.isColumnVector()) && (m_ncols != right.m_nrows)) {
        std::ostringstream oss;
        oss << "Matrix::operator* : (*this) (" << getNrows() << "x" << getNcols()
                << ") and RHS (" << right.getNrows() << "x" << right.getNcols()
                << ") do not have compatible dimensions";
        throw std::invalid_argument(oss.str().c_str());
    }
    Matrix result;
    switch (m_type) {
        case MATRIX_DENSE: // Left-hand side matrix is dense            
            result = multiplyLeftDense(right);
            break;
        case MATRIX_DIAGONAL:
            result = multiplyLeftDiagonal(right);
            break;
        case MATRIX_SYMMETRIC:
            result = multiplyLeftSymmetric(right);
            break;
        case MATRIX_SPARSE:
            result = multiplyLeftSparse(right);
            break;
        case MATRIX_LOWERTR:
            throw std::logic_error("Lower triangular multiplication not implemented yet");
        default:
            throw std::logic_error("unsupported");
    }
    return result;
}

Matrix & Matrix::operator=(const Matrix & right) {
    // Check for self-assignment!
    if (this == &right) {// Same object?
        return *this; // Yes, so skip assignment, and just return *this.
    }
    /* make sure shallow copies remain shallow */
    m_delete_data = (right.m_type != Matrix::MATRIX_SPARSE);
    m_ncols = right.m_ncols;
    m_nrows = right.m_nrows;
    m_type = right.m_type;
    m_triplet = NULL;
    m_sparse = NULL;
    m_dense = NULL;


    /* 
     * copy m_data only if 
     * (i)  the matrix is not sparse
     * (ii) the matrix is not shallow
     */
    m_dataLength = right.m_dataLength;
    m_data = NULL;
    if (right.m_type != MATRIX_SPARSE) {
        m_data = new double[m_dataLength];
        m_delete_data = true;
    }
    m_transpose = right.m_transpose;
    m_sparseStorageType = right.m_sparseStorageType;

    if (m_type == MATRIX_SPARSE) {
        if (right.m_triplet != NULL) {
            m_triplet = cholmod_copy_triplet(right.m_triplet, Matrix::cholmod_handle());
        }
        if (right.m_sparse != NULL) {
            m_sparse = cholmod_copy_sparse(right.m_sparse, Matrix::cholmod_handle());
        }
        if (right.m_dense != NULL) {
            m_dense = cholmod_copy_dense(right.m_dense, Matrix::cholmod_handle());
        }
    }
    if (Matrix::MATRIX_SPARSE != right.getType() && right.m_data != NULL) {
#ifdef USE_LIBS
        cblas_dcopy(m_dataLength, right.m_data, 1, m_data, 1);
#else
        for (int i = 0; i < m_dataLength; i++) {
            m_data[i] = right[i];
        }
#endif
    }
    return *this;
}

/********* PRIVATE METHODS ************/
void Matrix::domm(const Matrix &right, Matrix & result) const {
    // multiply with LHS being dense
    double t;
    for (size_t j = 0; j < right.m_ncols; j++) {
        for (size_t i = 0; i < m_nrows; i++) {
            t = 0.0;
            for (size_t k = 0; k < m_ncols; k++) {
                if (!(right.getType() == MATRIX_LOWERTR && k < j)) {
                    t += get(i, k) * right.get(k, j);
                }
            }
            result.set(i, j, t);
        }
    }
}

void Matrix::domm(Matrix& C, double alpha, Matrix& A, Matrix& B, double gamma) {
    // multiply with A being dense
    double t;
    for (size_t j = 0; j < B.m_ncols; j++) {
        for (size_t i = 0; i < C.m_nrows; i++) {
            t = 0.0;
            for (size_t k = 0; k < C.m_ncols; k++) {
                if (!(B.getType() == MATRIX_LOWERTR && k < j)) {
                    t += A.get(i, k) * B.get(k, j);
                }
            }
            C._addIJ(i, j, alpha*t, gamma);
        }
    }
}

Matrix Matrix::multiplyLeftDense(const Matrix & right) const {
    if (MATRIX_DENSE == right.m_type) { // RHS is also dense
        Matrix result(m_nrows, right.m_ncols);
#ifdef USE_LIBS
        cblas_dgemm(CblasColMajor,
                m_transpose ? CblasTrans : CblasNoTrans,
                right.m_transpose ? CblasTrans : CblasNoTrans,
                m_nrows, right.m_ncols, m_ncols, 1.0, m_data, m_transpose ? m_ncols : m_nrows,
                right.m_data, right.m_transpose ? right.m_ncols : right.m_nrows, 0.0,
                result.m_data, m_nrows);
#else
        domm(right, result);
#endif
        return result;
    } else if (MATRIX_DIAGONAL == right.m_type) { // {DENSE} * {DIAGONAL} = {DENSE} - RHS is diagonal
        Matrix result(getNrows(), getNcols());
        for (size_t j = 0; j < getNcols(); j++) {
            for (size_t i = 0; i < getNrows(); i++) {
                result.set(i, j, (!m_transpose ? m_data[i + j * m_nrows] : m_data[j + i * m_ncols]) * right.m_data[j]);
            }
        }
        return result;
    } else if (MATRIX_SYMMETRIC == right.m_type || MATRIX_LOWERTR == right.m_type) {
        Matrix result(m_nrows, right.m_ncols, Matrix::MATRIX_DENSE);
        domm(right, result);
        return result;
    } else { /* {DENSE} * {SPARSE} =  */
        /*
         * Trick: For D: dense, S: sparse it is
         * D * S = (S' * D')'
         */
        Matrix tempSparse(right);
        (const_cast<Matrix&> (*this)).transpose(); /*  D'  */
        tempSparse.transpose(); /*  S'  */
        Matrix r = tempSparse * (const_cast<Matrix&> (*this)); /*  r = S' * D'   */
        r.transpose(); /*  r'  */
        (const_cast<Matrix&> (*this)).transpose(); /*  D  */
        return r;
    }
}

Matrix Matrix::multiplyLeftSymmetric(const Matrix & right) const {
    // multiply when the LHS is symmetric    
    Matrix result(m_nrows, right.m_ncols);
    if (right.isColumnVector()) {
#ifdef USE_LIBS
        cblas_dspmv(CblasColMajor,
                CblasLower,
                m_nrows, 1.0, m_data,
                right.m_data, 1,
                0.0, result.m_data, 1);
        return result;
#endif
    }
    domm(right, result);
    return result;
}

Matrix Matrix::multiplyLeftDiagonal(const Matrix & right) const {
    // multiply when the LHS is diagonal
    Matrix result(m_nrows, right.m_ncols, right.m_type);
    for (size_t i = 0; i < m_nrows; i++) {
        if (MATRIX_SYMMETRIC == right.m_type) {
            for (size_t j = i; j < right.m_ncols; j++) {
                result.set(i, j, m_data[i] * right.get(i, j));
            }
        } else if (MATRIX_DENSE == right.m_type) {
            for (size_t j = 0; j < right.m_ncols; j++) {
                result.set(i, j, m_data[i] * right.get(i, j));
            }
        } else if (MATRIX_DIAGONAL == right.m_type) {
            result.set(i, i, m_data[i] * right.m_data[i]);
        } else if (MATRIX_LOWERTR == right.m_type) {
            for (size_t j = 0; j <= i; j++) {

                result.set(i, j, m_data[i] * right.get(i, j));
            }
        }
    }
    return result;
}

Matrix Matrix::multiplyLeftSparse(Matrix & right) {
    if (right.m_type == MATRIX_SPARSE) {
        // RHS is sparse
        if (m_sparse == NULL) {
            _createSparse();
        }
        if (right.m_sparse == NULL) {
            right._createSparse();
        }
        cholmod_sparse *r;
        r = cholmod_ssmult(
                (isColumnVector() && right.isColumnVector()) ? cholmod_transpose(m_sparse, 1, Matrix::cholmod_handle()) : m_sparse,
                right.m_sparse,
                0,
                true,
                false,
                Matrix::cholmod_handle());
        Matrix result(true);
        if (isColumnVector() && right.isColumnVector()) { /* Sparse-sparse dot product */
            result = Matrix(1, 1, Matrix::MATRIX_SPARSE);
        } else {
            result = Matrix(m_nrows, right.m_ncols, Matrix::MATRIX_SPARSE);
        }
        result.m_sparse = r;
        result._createTriplet();
        return result;
    } else if (right.m_type == MATRIX_DENSE) { /* SPRASE * DENSE */
        // RHS is dense
        Matrix result(getNrows(), right.getNcols());

        if (m_triplet != NULL && m_sparse == NULL)
            _createSparse();

        double alpha[2] = {1.0, 0.0};
        double beta[2] = {0.0, 0.0};

        if (right.m_dense == NULL) { /* Prepare right.m_dense */
            right.m_dense = cholmod_allocate_dense(
                    right.getNrows(),
                    right.getNcols(),
                    right.getNrows(),
                    CHOLMOD_REAL,
                    Matrix::cholmod_handle());

            if (!right.m_transpose) {
                right.m_dense->x = right.m_data;
            } else {
                /* Store RHS as transpose into right.m_dense */
                for (size_t i = 0; i < right.getNrows(); i++) { /* SPARSE x DENSE' */
                    for (size_t j = 0; j < right.getNcols(); j++) {
                        (static_cast<double*> (right.m_dense->x))[i + j * right.getNrows()] =
                                right.get(i, j);
                    }
                }
            }
        }

        bool dotProd = isColumnVector() && right.isColumnVector();
        result.m_dense = cholmod_allocate_dense(
                dotProd ? 1 : result.getNrows(),
                dotProd ? 1 : result.getNcols(),
                dotProd ? 1 : result.getNrows(),
                CHOLMOD_REAL,
                Matrix::cholmod_handle());
        if (m_transpose) {
            this->transpose();
        }
        cholmod_sdmult(
                m_sparse,
                dotProd ? true : m_transpose,
                alpha,
                beta,
                right.m_dense,
                result.m_dense,
                Matrix::cholmod_handle()
                );
        if (m_transpose) {
            this->transpose();
        }
        for (size_t k = 0; k < result.length(); k++) {
            result.m_data[k] = (static_cast<double*> (result.m_dense->x))[k];
        }
        return result;
    } else if (right.m_type == MATRIX_DIAGONAL) { // SPARSE * DIAGONAL = SPARSE
        Matrix result(*this); // COPY [result := right]
        for (size_t k = 0; k < result.m_triplet->nnz; k++) {
            int j_ = (static_cast<int*> (result.m_triplet->j))[k];
            (static_cast<double*> (result.m_triplet->x))[k] *= right.m_data[j_];
        }
        return result;
    } else {
        //LCOV_EXCL_START
        throw std::invalid_argument("SPARSE * {SYMMETRIC/LOWER/UPPER TRIANGUAL}: not supported");
        //LCOV_EXCL_STOP
    }
}

bool Matrix::indexWithinBounds(size_t i, size_t j) const {

    return (i < m_nrows && j < m_ncols) && !(m_type == MATRIX_LOWERTR && i < j);
}

Matrix::MatrixType Matrix::getType() const {

    return m_type;
}

void Matrix::init(size_t nr, size_t nc, MatrixType mType) {
    this -> m_transpose = false;
    this -> m_ncols = nc;
    this -> m_nrows = nr;
    this -> m_type = mType;
    this -> m_data = NULL;
    this -> m_delete_data = true;
    this -> m_triplet = NULL;
    this -> m_sparse = NULL;
    this -> m_dense = NULL;
    switch (m_type) {
        case MATRIX_DENSE:
            m_dataLength = nc * nr;
            m_data = new double[m_dataLength]();
            break;
        case MATRIX_DIAGONAL:
            if (nc != nr) {
                //LCOV_EXCL_START
                throw std::invalid_argument("Diagonal matrices must be square!!!");
                //LCOV_EXCL_STOP
            }
            m_dataLength = nc;
            m_data = new double[m_dataLength]();
            break;
        case MATRIX_LOWERTR:
        case MATRIX_SYMMETRIC:
            if (nc != nr) {
                //LCOV_EXCL_START
                throw std::invalid_argument("Lower triangular and symmetric matrices must be square!!!");
                //LCOV_EXCL_STOP
            }
            m_dataLength = nc * (nc + 1) / 2;
            m_data = new double[m_dataLength]();
            break;
        case MATRIX_SPARSE:
            m_data = NULL;
            break;
        default:
            throw std::logic_error("unsupported");
    }

}

void Matrix::_createSparse() {
    if (m_triplet != NULL) { // from triplets
        m_sparse = cholmod_triplet_to_sparse(m_triplet, m_triplet->nzmax, Matrix::cholmod_handle());
        return;
    }
    if (m_dense != NULL) { // from dense
        m_sparse = cholmod_dense_to_sparse(m_dense, true, Matrix::cholmod_handle());
        return;
    }
}

void Matrix::_createTriplet() {
    _createSparse();
    if (m_sparse != NULL) { /* make triplets from sparse */
        m_triplet = cholmod_sparse_to_triplet(m_sparse, Matrix::cholmod_handle());
    }
}

bool Matrix::isSymmetric() const {
    return (m_nrows == m_ncols) && ((Matrix::MATRIX_SYMMETRIC == m_type)
            || (Matrix::MATRIX_SPARSE == m_type && m_triplet != NULL && m_triplet->stype != 0)
            || (Matrix::MATRIX_DIAGONAL == m_type));
}

Matrix& operator*=(Matrix& obj, double alpha) {
    if (obj.m_type != Matrix::MATRIX_SPARSE) {
        assert(obj.m_data != NULL);
        cblas_dscal(obj.m_dataLength, alpha, obj.m_data, 1);
    } else {
        obj._createTriplet();
        obj.m_sparse = NULL;
        obj.m_dense = NULL;
        for (size_t k = 0; k < obj.m_triplet->nnz; k++) {
            (static_cast<double*> (obj.m_triplet->x))[k] *= alpha;
        }
    }
    return obj;
}

Matrix operator*(double alpha, Matrix& obj) {
    Matrix M(obj);
    M *= alpha;
    return (M);
}

std::string Matrix::getTypeString() const {
    //LCOV_EXCL_START
    const char names[5][10] = {
        "dense",
        "sparse",
        "diagonal",
        "lower",
        "symmetric"
    };

    int i = static_cast<int> (getType());
    std::string s = std::string(names[i]);
    return s;
    //LCOV_EXCL_STOP
}

Matrix Matrix::submatrixCopy(size_t row_start, size_t row_end, size_t col_start, size_t col_end) {
    /* DONE: Fully tested */
    //LCOV_EXCL_START
    if (row_end < row_start || col_end < col_start) {
        throw std::out_of_range("Matrix::submatrixCopy:: start > end is not allowed");
    }
    if (row_end > m_nrows) {
        throw std::out_of_range("Matrix::submatrixCopy:: row_end > total number of rows");
    }
    if (col_end > m_ncols) {
        throw std::out_of_range("Matrix::submatrixCopy:: col_end > total number of columns");
    }
    //LCOV_EXCL_STOP
    size_t rows = row_end - row_start + 1;
    size_t cols = col_end - col_start + 1;

    Matrix M(rows, cols, m_type);
    if (m_type == Matrix::MATRIX_DENSE) {
        /*
         * void dlacpy_( 
         *      char* uplo, 
         *      lapack_int* m, 
         *      lapack_int* n, 
         *      const double* a,
         *      lapack_int* lda, 
         *      double* b, 
         *      lapack_int* ldb );
         * 
         */
        dlacpy_(const_cast<char*> ("A"),
                reinterpret_cast<int*> (m_transpose ? &cols : &rows),
                reinterpret_cast<int*> (m_transpose ? &rows : &cols),
                m_data + (m_transpose ? row_start * m_ncols + col_start : row_start + col_start * m_nrows),
                const_cast<int*> (reinterpret_cast<const int*> (m_transpose ? &m_ncols : &m_nrows)),
                M.m_data,
                reinterpret_cast<int*> (m_transpose ? &cols : &rows));
        M.m_transpose = m_transpose;
    } else if (m_type == Matrix::MATRIX_SPARSE) {
        int * rs = new int[rows];
        int * cs = new int[cols];
        for (size_t i = 0; i <= rows; i++) {
            rs[i] = row_start + i;
        }
        for (size_t j = 0; j <= cols; j++) {
            cs[j] = col_start + j;
        }
        this->_createSparse();
        // int nnz = std::max(1, (int) (rows * cols / 20));
        cholmod_sparse * sp; // = cholmod_allocate_sparse(rows, cols, nnz, 1, 1, m_sparse->stype, m_sparse->xtype, cholmod_handle());
        sp = cholmod_submatrix(m_sparse, rs, rows, cs, cols, 1, 1, Matrix::cholmod_handle());
        M.m_sparse = sp;
        M._createTriplet();
    } else {
        //LCOV_EXCL_START
        throw std::logic_error("Matrix::submatrixCopy is available only for MATRIX_DENSE and MATRIX_SPARSE type matrices.");
        //LCOV_EXCL_STOP
    }
    return M;
}

Matrix Matrix::multiplySubmatrix(
        Matrix& right,
        const size_t left_row_start,
        const size_t left_row_end,
        const size_t left_col_start,
        const size_t left_col_end,
        const size_t right_row_start,
        const size_t right_row_end,
        const size_t right_col_start,
        const size_t right_col_end) {

    if (MATRIX_DENSE == m_type && MATRIX_DENSE == right.m_type) {
        size_t left_cols = left_col_end - left_col_start + 1;
        size_t left_rows = left_row_end - left_row_start + 1;

        size_t right_cols = right_col_end - right_col_start + 1;
        size_t right_rows = right_row_end - right_row_start + 1;

        if (left_cols != right_rows) {
            //LCOV_EXCL_START
            throw std::invalid_argument("Dimension of sub-matrix mismatch (left_cols!=right_rows)");
            //LCOV_EXCL_STOP
        }

        size_t left_start_idx =
                m_transpose
                ? left_row_start * m_ncols + left_col_start
                : left_row_start + left_col_start * m_nrows;
        size_t right_start_idx =
                right.m_transpose
                ? right_row_start * right.m_ncols + right_col_start
                : right_row_start + right_col_start * right.m_nrows;

        Matrix result(left_rows, right_cols, MATRIX_DENSE);

        cblas_dgemm(
                CblasColMajor,
                m_transpose ? CblasTrans : CblasNoTrans,
                right.m_transpose ? CblasTrans : CblasNoTrans,
                left_rows,
                right_cols,
                left_cols,
                1.0,
                m_data + left_start_idx,
                m_transpose ? m_ncols : m_nrows,
                right.m_data + right_start_idx,
                right.m_transpose ? right.m_ncols : right.m_nrows,
                0.0,
                result.m_data,
                left_rows);
        return result;

    } else {
        //LCOV_EXCL_START
        throw std::logic_error("Matrix::multiplySubmatrix is implemented only for dense matrix multiplication.");
        //LCOV_EXCL_STOP
    }
}

void Matrix::toggle_diagonal() {
    if (m_type != MATRIX_DENSE && m_type != MATRIX_DIAGONAL) { /* neither dense nor diagonal: unsupported. */
        throw std::invalid_argument("Only dense vectors and diagonal matrices are supported");
    }
    if (!isColumnVector() && m_type != MATRIX_DIAGONAL) {
        throw std::invalid_argument("Can only be applied to column vectors and diagonal matrices");
    }
    m_type = isColumnVector() ? MATRIX_DIAGONAL : MATRIX_DENSE;
    if (isColumnVector()) {
        m_ncols = m_nrows;
    } else {
        m_ncols = 1;
    }
}

Matrix::Matrix(bool shallow) {
    m_data = NULL;
    m_delete_data = false;
    m_nrows = 0;
    m_ncols = 0;
    m_dataLength = 0;
    m_type = MATRIX_DENSE;
    m_transpose = false;
    m_triplet = NULL;
    m_sparse = NULL;
    m_dense = NULL;
    m_sparseStorageType = CHOLMOD_TYPE_TRIPLET;
}

int Matrix::add(Matrix& C, double alpha, Matrix& A, double gamma) {
    // A and C must have compatible dimensions
    if (C.getNcols() != A.getNcols() || C.getNrows() != A.getNrows()) {
        throw std::invalid_argument("LHS and RHS do not have compatible dimensions");
    }
    // C := gamma * C + alpha * A
    int status;
    switch (C.getType()) {
        case MATRIX_DENSE: /* DENSE += ? */
            status = generic_add_helper_left_dense(C, alpha, A, gamma);
            break;
        case MATRIX_SYMMETRIC: /* SYMMETRIC += ? */
            status = generic_add_helper_left_symmetric(C, alpha, A, gamma);
            break;
        case MATRIX_LOWERTR: /* LOWER TRIANGULAR += ? */
            status = generic_add_helper_left_lower_tri(C, alpha, A, gamma);
            break;
        case MATRIX_DIAGONAL: /* DIAGONAL += ? */
            status = generic_add_helper_left_diagonal(C, alpha, A, gamma);
            break;
        case MATRIX_SPARSE: /* SPARSE += ? */
            status = generic_add_helper_left_sparse(C, alpha, A, gamma);
            break;
        default:
            status = ForBESUtils::STATUS_UNDEFINED_FUNCTION;
            break;
    }
    return status;
}

int Matrix::generic_add_helper_left_dense(Matrix& C, double alpha, Matrix& A, double gamma) {
    /*    HEREAFTER, C is DENSE!     */

    Matrix::MatrixType type_of_A = A.getType();

    bool is_gamma_one = (std::abs(gamma - 1.0) < std::numeric_limits<double>::epsilon());

    if (type_of_A == MATRIX_DENSE && C.m_transpose == A.m_transpose) {
        if (is_gamma_one) {
            cblas_daxpy(A.length(), alpha, A.m_data, 1, C.m_data, 1);
            return ForBESUtils::STATUS_OK;
        } else {
            for (size_t i = 0; i < A.length(); i++) {
                C.m_data[i] = (gamma * C.m_data[i]) + (alpha * A.m_data[i]);
            }
        }
    } else if (type_of_A == MATRIX_DIAGONAL) { /* DENSE + DIAGONAL */
        for (size_t i = 0; i < A.getNrows(); i++) {
            C._addIJ(i, i, alpha * A.m_data[i], gamma);
        }
    } else if (type_of_A == MATRIX_LOWERTR) { /* DENSE + LOWER */
        /* TODO This is to be tested!!! */
        for (size_t i = 0; i < A.getNrows(); i++) {
            if (!A.m_transpose) {
                for (size_t j = 0; j <= i; j++) {
                    C._addIJ(i, j, alpha * A.get(i, j), gamma);
                }
            } else {
                for (size_t j = i; j < A.m_ncols; j++) {
                    C._addIJ(i, j, alpha * A.get(i, j), gamma);
                }
            }
        }
    } else if (type_of_A == MATRIX_SPARSE) {
        A._createTriplet();
        assert(A.m_triplet != NULL);
        if (!is_gamma_one) {
            C *= gamma;
        }
        for (size_t k = 0; k < A.m_triplet->nnz; k++) {
            int i_ = (static_cast<int *> (A.m_triplet->i))[k];
            int j_ = (static_cast<int *> (A.m_triplet->j))[k];
            if (A.m_transpose) std::swap(i_, j_);
            C._addIJ(i_, j_, alpha * (static_cast<double *> (A.m_triplet->x))[k]);
        }
    } else { /* Symmetric and Dense+Dense' or Dense'+Dense (not of same transpose type) */
        for (size_t i = 0; i < A.getNrows(); i++) {
            for (size_t j = 0; j < A.getNcols(); j++) {
                C._addIJ(i, j, alpha * A.get(i, j), gamma);
            }
        }
    }
    return ForBESUtils::STATUS_OK;
}

int Matrix::generic_add_helper_left_symmetric(Matrix& C, double alpha, Matrix& A, double gamma) {
    int status = ForBESUtils::STATUS_OK;
    Matrix::MatrixType type_of_A = A.getType();
    size_t ncols = A.getNcols();
    size_t nrows = A.getNrows();
    bool is_gamma_one = (std::abs(gamma - 1.0) < std::numeric_limits<double>::epsilon());
    if (type_of_A == MATRIX_SYMMETRIC) {
        if (is_gamma_one) {
            cblas_daxpy(A.length(), alpha, A.m_data, 1, C.m_data, 1);
        } else {
            for (size_t i = 0; i < A.length(); i++) {
                C.m_data[i] = (gamma * C.m_data[i]) + (alpha * A.m_data[i]);
            }
        }
    } else if (type_of_A == MATRIX_DIAGONAL) { /* SYMMETRIC + DIAGONAL */
        for (size_t i = 0; i < ncols; i++) {
            size_t idx = i + C.m_nrows * i - i * (i + 1) / 2;
            double val = C.m_data[idx];
            C.m_data[idx] = gamma * val + alpha * A.m_data[i];
        }
    } else if (type_of_A == MATRIX_LOWERTR || type_of_A == MATRIX_DENSE) { /* SYMMETRIC + LOWER_TRI/DENSE = DENSE */
        C.m_dataLength = ncols * nrows; /* SYMMETRIC + DENSE = DENSE     */
        double * newData = new double[C.m_dataLength]; // realloc data
        for (size_t i = 0; i < nrows; i++) {
            for (size_t j = 0; j < ncols; j++) {
                newData[i + j * nrows] = gamma * C.get(i, j); // load data (recast into full storage format)
                if ((C.m_type == MATRIX_LOWERTR && type_of_A == MATRIX_DENSE) ? i >= j : true) {
                    newData[i + j * nrows] += alpha * A.get(i, j); // add RHS elements
                }
            }
        }
        double * temp_data;
        temp_data = static_cast<double *> (realloc(C.m_data, C.m_dataLength * sizeof (double))); // reallocate memory
        if (temp_data == NULL) { /* could not allocate using realloc */
            delete[] C.m_data;
            if (newData != NULL) {
                delete[] newData;
            }
            throw std::bad_alloc();
        } else if (temp_data == C.m_data) { /* nothing happened */
            temp_data = NULL;
        } else { /* realloc succeeded */
            C.m_data = temp_data; /* have m_data point to temp_data */
            temp_data = NULL; /* nullify temp_data */
        }
        C.m_data = static_cast<double *> (memcpy(C.m_data, newData, C.m_dataLength * sizeof (double))); // copy newData to m_data
        delete[] newData;
        C.m_type = MATRIX_DENSE;
        status = ForBESUtils::STATUS_HAD_TO_REALLOC;
    } else if (type_of_A == MATRIX_SPARSE) { /* SYMMETRIC + SPARSE */
        C.m_dataLength = ncols * nrows;
        double * newData = new double[C.m_dataLength];
        C.m_data = static_cast<double *> (realloc(C.m_data, C.m_dataLength * sizeof (double))); // reallocate memory

        for (size_t i = 0; i < nrows; i++) {
            for (size_t j = 0; j < ncols; j++) {
                newData[i + j * nrows] = C.get(i, j); /* restructure symmetric C data into dense */
            }
        }

        C.m_data = static_cast<double *> (memcpy(C.m_data, newData, C.m_dataLength * sizeof (double))); // copy newData to m_data

        A._createTriplet();
        assert(A.m_triplet != NULL);

        C.m_type = MATRIX_DENSE;

        for (size_t k = 0; k < A.m_triplet->nnz; k++) {
            int i_ = (static_cast<int*> (A.m_triplet->i))[k];
            int j_ = (static_cast<int*> (A.m_triplet->j))[k];
            C._addIJ(A.m_transpose ? j_ : i_, A.m_transpose ? i_ : j_, alpha * (static_cast<double*> (A.m_triplet->x))[k], gamma);
        }

        delete[] newData;
        status = ForBESUtils::STATUS_HAD_TO_REALLOC;
    }
    return status;
}

int Matrix::generic_add_helper_left_sparse(Matrix& C, double alpha, Matrix& A, double gamma) {
    int status = ForBESUtils::STATUS_OK;
    Matrix::MatrixType type_of_A = A.getType();
    size_t ncols = A.getNcols();
    size_t nrows = A.getNrows();

    if (type_of_A == MATRIX_SPARSE) {
        double __gamma_t[1] = {gamma};
        double __alpha_t[1] = {alpha};


        /*
         * Comments on the following two blocks of code:
         * IMPORTANT: READ ME
         * Don't forget that m_triplet is not really transposed - see #transpose().
         * When, at this point, m_sparse is NULL, the sparse data are stored in 
         * m_triplet. 
         * 
         * Recall that the code in #transpose() regarding sparse matrices is merely:
         * 
         * if (m_type == MATRIX_SPARSE) {
         *      _createSparse();
         *      m_sparse = cholmod_transpose(m_sparse, 1, cholmod_handle());
         * }
         * 
         * It doesn't modify m_triplet whatsoever. Now, take a look at #get():
         * 
         * double val = 0.0;
         * int i_ = m_transpose ? j : i;
         * int j_ = m_transpose ? i : j;
         * for (size_t k = 0; k < m_triplet->nnz; k++) {
         *     if (i_ == (static_cast<int*> (m_triplet->i))[k]) {
         *         if (j_ == (static_cast<int*> (m_triplet->j))[k]) {
         *             val = (static_cast<double*> (m_triplet->x))[k];
         *             break;
         *         }
         *     }
         * }
         * 
         * We swap ::i and ::j to get the correct entry and this is why the output
         * stream operator works as well (it uses #get(size_t, size_t)).
         * 
         * ERGO: If we create m_sparse out of the non-transposed m_triplet, hell
         * will brake loose. 
         */

        if (C.m_sparse == NULL) {
            C._createSparse(); /* CHOLMOD_SPARSE: create it 
                                from other representations */
            if (C.m_transpose) {
                C.m_sparse = cholmod_transpose(C.m_sparse, 1, cholmod_handle());
            }
        }

        if (A.m_sparse == NULL) {
            A._createSparse(); /* Likewise for the RHS */
            if (A.m_transpose) {
                A.m_sparse = cholmod_transpose(A.m_sparse, 1, cholmod_handle());
            }
        }

        /*
         * Don't store C as transpose any more - it will not have all values 
         * in the right order.
         */
        if (C.m_transpose) {
            C.m_transpose = false;
        }

        C.m_sparse = cholmod_add(
                C.m_sparse,
                A.m_sparse,
                __gamma_t,
                __alpha_t,
                true,
                true,
                Matrix::cholmod_handle()); /* Use cholmod_add to compute the sum C := gamma * C + alpha * A */

        C.m_triplet = cholmod_sparse_to_triplet(
                C.m_sparse,
                Matrix::cholmod_handle()); /* Update the triplet of the result (optional) */


    } else if (type_of_A == MATRIX_DIAGONAL) { /* SPARSE + DIAGONAL */
        C._createTriplet();
        if (C.m_triplet != NULL) {
            for (size_t i = 0; i < nrows; i++) {
                C._addIJ(i, i, alpha * A.get(i, i), gamma);
            }
        }
    } else if (type_of_A == MATRIX_LOWERTR) { /* SPARSE + LOWER TRI = SPARSE */
        C._createTriplet();
        if (C.m_triplet != NULL) {
            for (size_t i = 0; i < nrows; i++) {
                for (size_t j = (A.m_transpose ? i : 0);
                        j <= (A.m_transpose ? A.m_ncols : i); j++) {
                    C._addIJ(i, j, A.get(i, j), gamma);
                }
            }
        }
    } else if (type_of_A == MATRIX_DENSE) { /* SPARSE + DENSE */
        C.m_type = MATRIX_DENSE; /* Sparse + Dense = Dense */
        C.m_dataLength = A.getNcols() * A.getNrows(); /* Space to be allocated for the dense result */
        C.m_data = new double[C.m_dataLength]; /* allocate space */
        for (size_t i = 0; i < C.getNrows(); i++) {
            for (size_t j = 0; j < C.getNcols(); j++) {
                C.set(i, j, gamma * A.get(i, j)); /* store the rhs on m_data (accounting for transpose) */
            }
        }
        //_createTriplet();
        if (C.m_triplet != NULL) {
            /* Add triplets */
            for (size_t k = 0; k < C.m_triplet->nnz; k++) {
                int i = (static_cast<int*> (C.m_triplet->i))[k];
                int j = (static_cast<int*> (C.m_triplet->j))[k];
                if (C.m_transpose) {
                    std::swap(i, j);
                }
                double v = (static_cast<double*> (C.m_triplet->x))[k];
                C._addIJ(i, j, v, alpha);
            }
        }
        status = ForBESUtils::STATUS_HAD_TO_REALLOC;
    } else if (type_of_A == MATRIX_SYMMETRIC) { /* SPARSE + SYMMETRIC (result is dense) */
        C.m_dataLength = C.m_nrows * C.m_ncols;
        C.m_data = new double[C.m_dataLength](); // reallocate memory
        C.m_type = MATRIX_DENSE;
        for (size_t i = 0; i < nrows; i++) {
            for (size_t j = 0; j < ncols; j++) {
                C.set(i, j, gamma * A.get(i, j)); // load symmetric data                
            }
        }
        C._createTriplet();
        for (size_t k = 0; k < C.m_triplet->nnz; k++) {
            C._addIJ((static_cast<int*> (C.m_triplet->i))[k],
                    (static_cast<int*> (C.m_triplet->j))[k],
                    (static_cast<double*> (C.m_triplet->x))[k],
                    alpha);
        }
        status = ForBESUtils::STATUS_HAD_TO_REALLOC;
    }
    return status;
}

int Matrix::generic_add_helper_left_diagonal(Matrix& C, double alpha, Matrix& A, double gamma) {
    bool is_gamma_one = (std::abs(gamma - 1.0) < std::numeric_limits<double>::epsilon());
    if (MATRIX_DIAGONAL == A.m_type) { /* Diagonal += Diagonal */
        if (is_gamma_one) {
            cblas_daxpy(A.length(), alpha, A.m_data, 1, C.m_data, 1);
        } else {
            for (size_t i = 0; i < A.length(); i++) {
                C.m_data[i] = (gamma * C.m_data[i]) + (alpha * A.m_data[i]);
            }
        }
    } else {
        throw std::logic_error("Diagonal + Non-diagonal: not supported yet!");
    }
    return ForBESUtils::STATUS_OK;
}

int Matrix::generic_add_helper_left_lower_tri(Matrix& C, double alpha, Matrix& A, double gamma) {
    bool is_gamma_one = (std::abs(gamma - 1.0) < std::numeric_limits<double>::epsilon());
    if (A.m_type == Matrix::MATRIX_LOWERTR) { /* LOWER + LOWER = LOWER */
        if (is_gamma_one) {
            cblas_daxpy(A.length(), alpha, A.m_data, 1, C.m_data, 1);
        } else {
            for (size_t i = 0; i < A.length(); i++) {
                C.m_data[i] = (gamma * C.m_data[i]) + (alpha * A.m_data[i]);
            }
        }
    } else if (A.m_type == Matrix::MATRIX_DIAGONAL) {
        for (size_t i = 0; i < A.m_nrows; i++) {
            C._addIJ(i, i, alpha * A.get(i, i), gamma);
        }
    } else {
        throw std::logic_error("Lower-triangular + Non-lower triangular): not supported yet!");
    }
    return ForBESUtils::STATUS_OK;
}

int Matrix::mult(Matrix& C, double alpha, Matrix& A, Matrix& B, double gamma) {
    // A and C must have compatible dimensions
    if (A.getNcols() != B.getNrows()) {
        std::ostringstream oss;
        oss << "A (" << A.getNrows() << "x" << A.getNcols()
                << ") and B (" << B.getNrows() << "x" << B.getNcols()
                << ") do not have compatible dimensions";
        throw std::invalid_argument(oss.str().c_str());
    }
    /* C must have proper dimensions */
    if (C.getNrows() != A.getNrows() || C.getNcols() != B.getNcols()) {
        std::ostringstream oss;
        oss << "C is " << C.getNrows() << "x" << C.getNrows()
                << ", but it should be " << A.getNrows() << "x"
                << B.getNcols();
        throw std::invalid_argument(oss.str().c_str());
    }
    // C := gamma * C + alpha * A * B
    int status = ForBESUtils::STATUS_UNDEFINED_FUNCTION;
    switch (A.getType()) {
        case MATRIX_DENSE: /* DENSE += ? */
            status = multiply_helper_left_dense(C, alpha, A, B, gamma);
            break;
        case MATRIX_SYMMETRIC: /* SYMMETRIC += ? */
            status = multiply_helper_left_symmetric(C, alpha, A, B, gamma);
            break;
        case MATRIX_LOWERTR: /* LOWER TRIANGULAR += ? */
            //            status = generic_add_helper_left_lower_tri(C, alpha, A, gamma);
            break;
        case MATRIX_DIAGONAL: /* DIAGONAL += ? */
            status = multiply_helper_left_diagonal(C, alpha, A, B, gamma);
            break;
        case MATRIX_SPARSE: /* SPARSE += ? */
            status = multiply_helper_left_sparse(C, alpha, A, B, gamma);
            break;
        default:
            break;
    }
    return status;
}

int Matrix::multiply_helper_left_dense(Matrix& C, double alpha, Matrix& A, Matrix& B, double gamma) {
    /* A is dense */
    int status = ForBESUtils::STATUS_OK;
    if (C.m_dataLength < C.getNrows() * C.getNcols()) {
        C = Matrix(C.getNrows(), C.getNcols(), Matrix::MATRIX_DENSE);
        status = ForBESUtils::STATUS_HAD_TO_REALLOC;
    }
    C.m_type = Matrix::MATRIX_DENSE;
    if (MATRIX_DENSE == B.m_type) { // B is also dense    
        cblas_dgemm(CblasColMajor,
                A.m_transpose ? CblasTrans : CblasNoTrans,
                B.m_transpose ? CblasTrans : CblasNoTrans,
                A.m_nrows,
                B.m_ncols,
                A.m_ncols,
                alpha,
                A.m_data,
                A.m_transpose ? A.m_ncols : A.m_nrows,
                B.m_data,
                B.m_transpose ? B.m_ncols : B.m_nrows,
                gamma,
                C.m_data,
                C.m_nrows);
        status = std::max(status, ForBESUtils::STATUS_OK);
    } else if (MATRIX_DIAGONAL == B.m_type) { // {DENSE} * {DIAGONAL} = {DENSE} - B is diagonal
        for (size_t j = 0; j < C.getNcols(); j++) {
            for (size_t i = 0; i < C.getNrows(); i++) {
                C._addIJ(i, j, alpha * A.get(i, j) * B.get(j, j), gamma);
            }
        }
        status = ForBESUtils::STATUS_OK;
    } else if (MATRIX_SYMMETRIC == B.m_type || MATRIX_LOWERTR == B.m_type) {
        domm(C, alpha, A, B, gamma);
        status = ForBESUtils::STATUS_OK;
    } else { /* {DENSE} * {SPARSE} =  */
        /*
         * Trick: For D: dense, S: sparse it is
         * D * S = (S' * D')'
         */
        //        Matrix tempSparse(right);
        //        (const_cast<Matrix&> (*this)).transpose(); /*  D'  */
        //        tempSparse.transpose(); /*  S'  */
        //        Matrix r = tempSparse * (const_cast<Matrix&> (*this)); /*  r = S' * D'   */
        //        r.transpose(); /*  r'  */
        //        (const_cast<Matrix&> (*this)).transpose(); /*  D  */
        //        return r;
    }
    return status;
}

int Matrix::multiply_helper_left_sparse(Matrix& C, double alpha, Matrix& A, Matrix& B, double gamma) {
    bool is_alpha_one = (std::abs(alpha - 1.0) < std::numeric_limits<double>::epsilon());
    bool is_gamma_zero = (std::abs(gamma) < std::numeric_limits<double>::epsilon());
    bool is_gamma_one = (std::abs(gamma) < std::numeric_limits<double>::epsilon());
    int status = ForBESUtils::STATUS_UNDEFINED_FUNCTION;
    if (B.m_type == MATRIX_SPARSE) {
        // RHS is sparse
        if (A.m_sparse == NULL) {
            A._createSparse();
        }
        if (B.m_sparse == NULL) {
            B._createSparse();
        }
        if (C.m_sparse == NULL) {
            C._createSparse();
        }
        cholmod_sparse *r; // r will store A * B
        r = cholmod_ssmult(
                (A.isColumnVector() && B.isColumnVector())
                ? cholmod_transpose(A.m_sparse, 1, Matrix::cholmod_handle())
                : A.m_sparse,
                B.m_sparse,
                0,
                true,
                false,
                Matrix::cholmod_handle()); // r = A*B

        /*
         * SCALE: r *= alpha (unless alpha == 1)
         */
        if (!is_alpha_one) {
            for (size_t j = 0; j < C.m_ncols; j++) {
                int p = (static_cast<int*> (r->p))[j];
                int pend = (r->packed == 1)
                        ? ((static_cast<int*> (r->p))[j + 1])
                        : p + (static_cast<int*> (r->nz))[j];
                for (; p < pend; p++) {
                    (static_cast<double*> (r->x))[p] *= alpha;
                }
            }
        }


        if (is_gamma_zero) {
            cholmod_triplet * r_to_triplet = cholmod_sparse_to_triplet(r, Matrix::cholmod_handle());
            // C := alpha * A * B = r
            C.m_sparse = r;
            C.m_triplet = r_to_triplet;
            status = ForBESUtils::STATUS_OK;
        } else {
            Matrix temp_r = Matrix(true);
            temp_r.m_nrows = C.getNrows();
            temp_r.m_ncols = C.getNcols();
            temp_r.m_sparse = cholmod_copy_sparse(r, cholmod_handle());
            temp_r.m_type = MATRIX_SPARSE;
            temp_r.m_triplet = cholmod_sparse_to_triplet(temp_r.m_sparse, Matrix::cholmod_handle());
            add(C, 1.0, temp_r, gamma);
            status = ForBESUtils::STATUS_OK;
        }
    } else if (B.m_type == MATRIX_DENSE) { /* C = gamma * C + alpha * SPRASE * DENSE */


        // RHS is dense
        //                Matrix result(A.getNrows(), B..getNcols());

        //        if (A.m_triplet != NULL && A.m_sparse == NULL)
        //            A._createSparse();
        //
        //        double alpha[2] = {1.0, 0.0};
        //        double beta[2] = {0.0, 0.0};
        //
        //        if (B.m_dense == NULL) { /* Prepare right.m_dense */
        //            B.m_dense = cholmod_allocate_dense(
        //                    B.getNrows(),
        //                    B.getNcols(),
        //                    B.getNrows(),
        //                    CHOLMOD_REAL,
        //                    Matrix::cholmod_handle());
        //
        //            if (!B.m_transpose) {
        //                B.m_dense->x = B.m_data;
        //            } else {
        //                /* Store RHS as transpose into right.m_dense */
        //                for (size_t i = 0; i < B.getNrows(); i++) { /* SPARSE x DENSE' */
        //                    for (size_t j = 0; j < B.getNcols(); j++) {
        //                        (static_cast<double*> (B.m_dense->x))[i + j * B.getNrows()] = B.get(i, j);
        //                    }
        //                }
        //            }
        //        }
        //
        //        bool dotProd = isColumnVector() && B.isColumnVector();
        //        C.m_dense = cholmod_allocate_dense(
        //                dotProd ? 1 : B.getNrows(),
        //                dotProd ? 1 : B.getNcols(),
        //                dotProd ? 1 : B.getNrows(),
        //                CHOLMOD_REAL,
        //                Matrix::cholmod_handle());

        //        if (m_transpose) {
        //            this->transpose();
        //        }
        //        cholmod_sdmult(
        //                m_sparse,
        //                dotProd ? true : m_transpose,
        //                alpha,
        //                beta,
        //                right.m_dense,
        //                result.m_dense,
        //                Matrix::cholmod_handle()
        //                );
        //        if (m_transpose) {
        //            this->transpose();
        //        }
        //        for (size_t k = 0; k < result.length(); k++) {
        //            result.m_data[k] = (static_cast<double*> (result.m_dense->x))[k];
        //        }


    } else if (B.m_type == MATRIX_DIAGONAL) { // += alpha * SPARSE * DIAGONAL
        Matrix A_temp(A); //  Compute A_temp = A * alpha;
        for (size_t k = 0; k < A.m_triplet->nnz; k++) {
            int j_ = (static_cast<int*> (A_temp.m_triplet->j))[k];
            static_cast<double*> (A_temp.m_triplet->x)[k] *= (alpha * B.get(j_, j_));
        }
        status = add(C, 1.0, A_temp, gamma);
    } else {
        //LCOV_EXCL_START
        throw std::invalid_argument("SPARSE * {SYMMETRIC/LOWER/UPPER TRIANGUAL}: not supported");
        //LCOV_EXCL_STOP
    }
    return status;
}

int Matrix::multiply_helper_left_diagonal(Matrix& C, double alpha, Matrix& A, Matrix& B, double gamma) {
    for (size_t i = 0; i < C.m_nrows; i++) {
        if (MATRIX_SYMMETRIC == B.m_type) {
            for (size_t j = i; j < B.m_ncols; j++) {
                C._addIJ(i, j, alpha * A.m_data[i] * B.get(i, j), gamma);
            }
        } else if (MATRIX_DENSE == B.m_type) {
            for (size_t j = 0; j < B.m_ncols; j++) {
                C._addIJ(i, j, alpha * A.m_data[i] * B.get(i, j), gamma);
            }
        } else if (MATRIX_DIAGONAL == B.m_type) {
            C._addIJ(i, i, alpha * A.m_data[i] * B.m_data[i], gamma);
        } else if (MATRIX_LOWERTR == B.m_type) {
            for (size_t j = 0; j <= i; j++) {
                C._addIJ(i, j, alpha * A.m_data[i] * B.get(i, j), gamma);
            }
        }
    }
    return ForBESUtils::STATUS_OK;
}

int Matrix::multiply_helper_left_symmetric(Matrix& C, double alpha, Matrix& A, Matrix& B, double gamma) {
    // multiply when the LHS is symmetric        
    if (B.isColumnVector()) {
        cblas_dspmv(CblasColMajor,
                CblasLower,
                A.m_nrows,
                alpha,
                A.m_data,
                B.m_data,
                1,
                gamma,
                C.m_data,
                1);
    } else {
        domm(C, alpha, A, B, gamma);
    }
    return ForBESUtils::STATUS_OK;
}
