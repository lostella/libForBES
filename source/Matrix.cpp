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
#include "MatrixFactory.h"
#include <assert.h>

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
    m_data = new double;
    *m_data = 0;
    m_type = MATRIX_DENSE;
    m_dataLength = 0;
    m_transpose = false;
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
    if (orig.m_type != MATRIX_SPARSE) {
        size_t n = orig.m_dataLength;
        if (n <= 0) {
            n = 1;
        }
        m_data = new double[n];
        for (size_t i = 0; i < n; i++) {
            m_data[i] = orig.m_data[i];
        }
        m_dataLength = orig.m_dataLength;
    }
    m_type = orig.m_type;
    if (m_type == MATRIX_SPARSE) {
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
}

/********* DENSTRUCTOR ************/
Matrix::~Matrix() {
    this -> m_ncols = 0;
    this -> m_nrows = 0;
    if (m_data != NULL) {
        delete[] m_data;
    }
    m_data = NULL;
    if (m_triplet != NULL) {
        cholmod_free_triplet(&m_triplet, Matrix::cholmod_handle());
        m_triplet = NULL;
    }
    if (m_sparse != NULL) {
        cholmod_free_sparse(&m_sparse, Matrix::cholmod_handle());
        m_sparse = NULL;
    }
    if (m_dense != NULL) {
        cholmod_free_dense(&m_dense, Matrix::cholmod_handle());
    }
}

/********* GETTERS/SETTERS ************/
size_t Matrix::getNcols() const {
    return m_ncols;
}

size_t Matrix::getNrows() const {
    return m_nrows;
}

double * const Matrix::getData() const {
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
    size_t new_size = nrows*ncols;
    if (new_size <= 0) {
        return -1;
    }
    if (new_size > length()) {
        return -2;
    }
    this -> m_nrows = nrows;
    this -> m_ncols = ncols;
    return 0;
}

double Matrix::get(const size_t i, const size_t j) const {
    if (isEmpty()) {
        throw std::out_of_range("Method get(size_t, size_t) applied to an empty matrix");
    }
    if (i < 0 || i >= getNrows() || j < 0 || j >= getNcols()) {
        throw std::out_of_range("Index out of range!");
    }
    if (m_type == MATRIX_DENSE) {
        return !m_transpose ? m_data[i + j * m_nrows] : m_data[j + i * m_ncols];
    } else if (m_type == MATRIX_DIAGONAL) {
        if (i == j) {
            return m_data[i];
        } else {
            return 0.0f;
        }
    } else if (m_type == MATRIX_SYMMETRIC) {
        size_t i_ = std::max(i, j);
        size_t j_ = std::min(i, j);
        return m_data[i_ + m_nrows * j_ - j_ * (j_ + 1) / 2];
    } else if (m_type == MATRIX_LOWERTR) {
        return m_transpose
                ? (j >= i) ? m_data[j + m_ncols * i - i * (i + 1) / 2] : 0.0f
                : (i >= j) ? m_data[i + m_nrows * j - j * (j + 1) / 2] : 0.0f;
    } else {
        /* if (m_type == MATRIX_SPARSE) */
        if (m_triplet == NULL) {
            throw std::logic_error("not supported yet");
        }
        double val = 0.0f;
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
    if (!indexWithinBounds(i, j)) {
        throw std::out_of_range("Index out of range!");
    }
    if (m_type == MATRIX_DENSE) {
        m_data[i + j * m_nrows] = v;
    } else if (m_type == MATRIX_DIAGONAL && i == j) {
        m_data[i] = v;
    } else if (m_type == MATRIX_LOWERTR) {
        m_data[i + m_nrows * j - j * (j + 1) / 2] = v;
    } else if (m_type == MATRIX_SYMMETRIC) {
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

        throw std::invalid_argument("Illegal operation");
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

double Matrix::quad(Matrix & x) {
    if (!x.isColumnVector()) {
        throw std::invalid_argument("Method `quadratic` can only be applied to vectors!");
    }
    if (getNcols() != getNrows()) {
        throw std::invalid_argument("Method `quadratic` can only be applied on square matrices!");
    }
    if (x.getNrows() != m_ncols) {
        throw std::invalid_argument("The argument of quad(Matrix&) is not of appropriate dimension.");
    }
    double result = 0.0;

    if (MATRIX_DENSE == m_type || MATRIX_LOWERTR == m_type) { /* DENSE or LOWER TRIANGULAR */
        for (size_t i = 0; i < m_nrows; i++) {
            for (size_t j = 0; j < m_ncols; j++) {
                result += x[i] * get(i, j) * x[j];
            }
        }
    } else if (MATRIX_DIAGONAL == m_type) { /* DIAGONAL */
        for (size_t i = 0; i < m_nrows; i++) {
            result += x[i] * x[i] * get(i, i);
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
            /* TODO: Implement quadFromSparse */
            throw std::logic_error("Quad on sparse matrix - no triplets found (not implemented yet)");
        }
    }
    result /= 2;

    return result;
}

double Matrix::quad(Matrix& x, Matrix & q) {
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
    double t = 0.0f;
    Matrix r;
    r = q*x;
    assert(!r.isEmpty());
    t = r.get(0, 0) + quad(x); /* q*x + (1/2) x'*Q*x */

    return t;
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
            int p = ((int*) obj.m_sparse->p)[j];
            int pend = (obj.m_sparse->packed == 1) ? (((int*) obj.m_sparse->p)[j + 1]) : p + ((int*) obj.m_sparse->nz)[j];
            for (; p < pend; p++) {
                os << "(" << ((int*) obj.m_sparse->i)[p] << "," << j << ")  : "
                        << std::setw(8) << std::setprecision(4) << ((double*) obj.m_sparse->x)[p] << std::endl;
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
    if (sub < 0 || sub >= length()) {
        throw std::out_of_range("Exception: Index out of range for Matrix");
    }
    return m_data[sub];
}

inline void Matrix::_addIJ(size_t i, size_t j, double a) {
    set(i, j, get(i, j) + a); // A(i,j) += a
}

inline void Matrix::_addD(Matrix& rhs) { /* DENSE += (?) */
    /* LHS is DENSE */
    assert(MATRIX_DENSE == m_type);
    if (rhs.m_type == MATRIX_DENSE && m_transpose == rhs.m_transpose) {
        /* RHS is also dense and of same transpose-type (D+D or D'+D') */
        vectorAdd(length(), m_data, rhs.m_data);
    } else if (rhs.getType() == MATRIX_DIAGONAL) { /* DENSE + DIAGONAL */
        for (size_t i = 0; i < m_nrows; i++) {
            _addIJ(i, i, rhs.get(i, i));
        }
    } else if (rhs.getType() == MATRIX_LOWERTR) { /* DENSE + LOWER */
        /* TODO This is to be tested!!! */
        for (size_t i = 0; i < m_nrows; i++) {
            if (!rhs.m_transpose) {
                for (size_t j = 0; j <= i; j++) {
                    _addIJ(i, j, rhs.get(i, j));
                }
            } else {
                for (size_t j = i; j < rhs.m_ncols; j++) {
                    _addIJ(i, j, rhs.get(i, j));
                }
            }
        }
    } else if (rhs.getType() == MATRIX_SPARSE) {
        rhs._createTriplet();
        assert(rhs.m_triplet != NULL);
        for (size_t k = 0; k < rhs.m_triplet->nnz; k++) {
            int i_ = ((int*) rhs.m_triplet->i)[k];
            int j_ = ((int*) rhs.m_triplet->j)[k];
            _addIJ(rhs.m_transpose ? j_ : i_, rhs.m_transpose ? i_ : j_, ((double*) rhs.m_triplet->x)[k]);
        }
    } else { /* Symmetric and Dense+Dense' or Dense'+Dense (not of same transpose type) */
        for (size_t i = 0; i < m_nrows; i++) {
            for (size_t j = 0; j < m_ncols; j++) {
                _addIJ(i, j, rhs.get(i, j));
            }
        }
    }
}

inline void Matrix::_addH(Matrix& rhs) { /* SYMMETRIC += (?) */
    assert(MATRIX_SYMMETRIC == m_type);
    if (rhs.m_type == MATRIX_SYMMETRIC) {
        vectorAdd(length(), m_data, rhs.m_data);
    } else if (rhs.m_type == MATRIX_DIAGONAL) { /* SYMMETRIC + DIAGONAL */
        for (size_t i = 0; i < m_ncols; i++) {
            _addIJ(i, i, rhs.get(i, i));
        }
    } else if (rhs.m_type == MATRIX_LOWERTR || rhs.m_type == MATRIX_DENSE) { /* SYMMETRIC + LOWER_TRI/DENSE = DENSE */
        m_dataLength = m_ncols * m_nrows; /* SYMMETRIC + DENSE = DENSE     */
        double * newData = new double[m_dataLength]; // realloc data
        for (size_t i = 0; i < m_nrows; i++) {
            for (size_t j = 0; j < m_ncols; j++) {
                newData[i + j * m_nrows] = get(i, j); // load data (recast into full storage format)
                if ((m_type == MATRIX_LOWERTR && rhs.m_type == MATRIX_DENSE) ? i >= j : true) {
                    newData[i + j * m_nrows] += rhs.get(i, j); // add RHS elements
                }
            }
        }
        m_data = (double*) realloc(m_data, m_dataLength * sizeof (double)); // reallocate memory
        m_data = (double*) memcpy(m_data, newData, m_dataLength * sizeof (double)); // copy newData to m_data
        delete newData;
        m_type = MATRIX_DENSE;
    } else if (rhs.m_type == MATRIX_SPARSE) { /* SYMMETRIC + SPARSE */
        m_dataLength = m_ncols * m_nrows;
        double * newData = new double[m_dataLength];
        m_data = (double *) realloc(m_data, m_dataLength * sizeof (double)); // reallocate memory

        for (size_t i = 0; i < m_nrows; i++) {
            for (size_t j = 0; j < m_ncols; j++) {
                newData[i + j * m_nrows] = get(i, j); /* restructure symmetric data into dense */
            }
        }

        m_data = (double *) memcpy(m_data, newData, m_dataLength * sizeof (double)); // copy newData to m_data

        rhs._createTriplet();
        assert(rhs.m_triplet != NULL);

        m_type = MATRIX_DENSE;

        for (size_t k = 0; k < rhs.m_triplet->nnz; k++) {
            int i_ = ((int*) rhs.m_triplet->i)[k];
            int j_ = ((int*) rhs.m_triplet->j)[k];
            _addIJ(rhs.m_transpose ? j_ : i_, rhs.m_transpose ? i_ : j_, ((double*) rhs.m_triplet->x)[k]);
        }

        delete newData;
    }
}

inline void Matrix::vectorAdd(size_t len, double* pV1, const double* pV2) {
#ifdef USE_LIBS 
    cblas_daxpy(len, 1.0, pV2, 1, pV1, 1); // data = data + right.data
#else
    for (int i = 0; i < len; i++)
        pV1[i] += pV2[i];
#endif
}

inline void Matrix::_addL(Matrix& rhs) { /* LOWER += (?) */
    assert(m_type == MATRIX_LOWERTR);
    if (rhs.m_type == Matrix::MATRIX_LOWERTR) { /* LOWER + LOWER = LOWER */
        vectorAdd(length(), m_data, rhs.m_data);
    } else if (rhs.m_type == Matrix::MATRIX_DIAGONAL) {
        for (size_t i = 0; i < m_nrows; i++) {
            _addIJ(i, i, rhs.get(i, i));
        }
    } else {
        throw std::logic_error("Lower-triangular + Non-lower triangular): not supported yet!");
    }
}

inline void Matrix::_addS(Matrix& rhs) { /* SPARSE += (?) */
    assert(m_type == MATRIX_SPARSE);
    if (rhs.m_type == MATRIX_SPARSE) {
        double alpha[1] = {1};

        if (m_sparse == NULL)
            _createSparse(); /* CHOLMOD_SPARSE: create it 
                                from other representations */

        if (rhs.m_sparse == NULL)
            rhs._createSparse(); /* Likewise for the RHS */


        m_sparse = cholmod_add(m_sparse,
                rhs.m_sparse,
                alpha,
                alpha,
                true,
                true,
                Matrix::cholmod_handle()); /* Use cholmod_add to compute the sum A+=B */

        m_triplet = cholmod_sparse_to_triplet(
                m_sparse,
                Matrix::cholmod_handle()); /* Update the triplet of the result (optional) */
    } else if (rhs.m_type == MATRIX_DIAGONAL) { /* SPARSE + DIAGONAL */
        _createTriplet();
        if (m_triplet != NULL) {
            for (size_t i = 0; i < getNrows(); i++) {
                _addIJ(i, i, rhs.get(i, i));
            }
        }
    } else if (rhs.m_type == MATRIX_LOWERTR) { /* SPARSE + LOWER TRI = SPARSE */
        _createTriplet();
        if (m_triplet != NULL) {
            for (size_t i = 0; i < m_nrows; i++) {
                for (size_t j = (rhs.m_transpose ? i : 0);
                        j <= (rhs.m_transpose ? rhs.m_ncols : i); j++) {
                    _addIJ(i, j, rhs.get(i, j));
                }
            }
        }
    } else if (rhs.m_type == MATRIX_DENSE) { /* For any type of non-sparse right-summands */
        /* Sparse + Dense = Dense */
        m_type = MATRIX_DENSE;
        m_dataLength = m_ncols * m_nrows;
        m_data = new double[m_ncols * m_nrows];
        for (size_t i = 0; i < m_nrows; i++) {
            for (size_t j = 0; j < m_ncols; j++) {
                set(i, j, rhs.get(i, j));
            }
        }
        _createTriplet();
        if (m_triplet != NULL) {
            /* Add triplets */
            int i = -1, j = -1;
            double v = 0.0;
            for (size_t k = 0; k < m_triplet->nnz; k++) {
                i = ((int*) m_triplet->i)[k];
                j = ((int*) m_triplet->j)[k];
                v = ((double*) m_triplet->x)[k];
                m_data[i + j * m_nrows] += v;
            }
        }
    } else if (rhs.m_type == MATRIX_SYMMETRIC) { /* SPARSE + SYMMETRIC (result is dense) */
        m_dataLength = m_nrows * m_ncols;
        m_data = new double[m_dataLength](); // reallocate memory
        m_type = MATRIX_DENSE;
        for (size_t i = 0; i < m_nrows; i++) {
            for (size_t j = 0; j < m_ncols; j++) {
                set(i, j, rhs.get(i, j)); // load symmetric data                
            }
        }
        _createTriplet();
        for (size_t k = 0; k < m_triplet->nnz; k++) {
            _addIJ(((int*) m_triplet->i)[k], ((int*) m_triplet->j)[k], ((double*) m_triplet->x)[k]);
        }

    }

}

void Matrix::_addX(Matrix& rhs) {
    assert(m_type == MATRIX_DIAGONAL);
    if (MATRIX_DIAGONAL == rhs.m_type) { /* Diagonal += Diagonal */
        vectorAdd(length(), m_data, rhs.m_data);
    } else {
        throw std::logic_error("Diagonal + Non-diagonal: not supported yet!");
    }
}

Matrix& Matrix::operator+=(Matrix & right) {
    if (m_ncols != right.m_ncols || m_nrows != right.m_nrows) {
        throw std::invalid_argument("Incompatible dimensions while using +=!");
    }

    if (&right == this) {
        *this *= 2.0;
        return *this;
    }

    switch (m_type) {
        case MATRIX_DENSE: /* DENSE += ? */
            _addD(right);
            break;
        case MATRIX_SYMMETRIC: /* SYMMETRIC += ? */
            _addH(right);
            break;
        case MATRIX_LOWERTR: /* LOWER TRIANGULAR += ? */
            _addL(right);
            break;
        case MATRIX_DIAGONAL: /* DIAGONAL += ? */
            _addX(right);
            break;
        case MATRIX_SPARSE: /* SPARSE += ? */
            _addS(right);
            break;
    }
    return *this;
}

Matrix & Matrix::operator-=(const Matrix & right) {
    if (m_ncols != right.m_ncols || m_nrows != right.m_nrows) {
        throw std::invalid_argument("Incompatible dimensions while using +=!");
    }
#ifdef USE_LIBS
    cblas_daxpy(length(), -1.0f, right.m_data, 1, m_data, 1); // data = data + right.data
#else
    for (int i = 0; i < length(); i++) {

        m_data[i] -= right[i];
    }
#endif
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

Matrix Matrix::operator-(const Matrix & right) const {
    if (this->getNrows() != right.getNrows() || this->getNcols() != right.getNcols()) {
        throw std::invalid_argument("Addition of matrices of incompatible dimensions!");
    }
    Matrix result(*this); // Make a copy of myself.      
    result -= right;

    return result;
}

Matrix Matrix::operator*(Matrix & right) {
    double t = 0.0f;
    if (!(getType() == Matrix::MATRIX_SPARSE && right.getType() == Matrix::MATRIX_SPARSE) &&
            isColumnVector() && right.isColumnVector() && length() == right.length()) {
        // multiplication of two column vectors = dot product
        Matrix r(1, 1);
        for (size_t i = 0; i < m_nrows * m_ncols; i++) {
            t += m_data[i] * right.m_data[i];
        }
        r[0] = t;
        return r;
    }
    if (!(isColumnVector() && right.isColumnVector()) && (m_ncols != right.m_nrows)) {
        throw std::invalid_argument("Incompatible dimensions!");
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
    }
    return result;
}

Matrix & Matrix::operator=(const Matrix & right) {
    // Check for self-assignment!
    if (this == &right) {// Same object?
        return *this; // Yes, so skip assignment, and just return *this.
    }
    m_ncols = right.m_ncols;
    m_nrows = right.m_nrows;
    m_type = right.m_type;
    if (right.m_type != MATRIX_SPARSE) {
        m_dataLength = right.m_dataLength;
        if (m_data != NULL) {
            delete m_data;
        }
        m_data = new double[m_dataLength];
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
    //std::cout << "multiplication - left = " << m_nrows << "x" << m_ncols << std::endl;
    for (size_t i = 0; i < m_nrows; i++) {
        for (size_t j = 0; j < right.m_ncols; j++) {
            t = 0.0;
            //std::cout << "For element C(" << i << "," << j << ")\n";
            for (size_t k = 0; k < m_ncols; k++) {
                //std::cout << "  += A(" << i << "," << k << ") * B(" << k << "," << i << ")" << std::endl;
                if (!(right.getType() == MATRIX_LOWERTR && k < j)) {

                    t += get(i, k) * right.get(k, j);
                }
            }
            result.set(i, j, t);
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
                m_nrows, right.m_ncols, m_ncols, 1.0f, m_data, m_transpose ? m_ncols : m_nrows,
                right.m_data, right.m_transpose ? right.m_ncols : right.m_nrows, 0.0f,
                result.m_data, m_nrows);
#else
        domm(right, result);
#endif
        return result;
    } else if (MATRIX_DIAGONAL == right.m_type) { // {DENSE} * {DIAGONAL} = {DENSE} - RHS is diagonal
        Matrix result(*this);
        for (size_t j = 0; j < m_ncols; j++) {
            for (size_t i = 0; i < m_nrows; i++) {
                result.set(i, j, result.get(i, j) * right.get(j, j));
            }
        }
        return result;
    } else if (MATRIX_SYMMETRIC == right.m_type || MATRIX_LOWERTR == right.m_type) {
        Matrix result(m_nrows, right.m_ncols, Matrix::MATRIX_DENSE);
        domm(right, result);
        return result;
    } else {

        throw std::logic_error("We apologize for not having implemented this functionality yet (Dense * Sparse)!");
    }
}

Matrix Matrix::multiplyLeftSymmetric(const Matrix & right) const {
    // multiply when the LHS is symmetric    
    Matrix result(m_nrows, right.m_ncols);
    if (right.isColumnVector()) {
#ifdef USE_LIBS
        cblas_dspmv(CblasColMajor, CblasLower,
                m_nrows, 1.0f, m_data,
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
                result.set(i, j, get(i, i) * right.get(i, j));
            }
        } else if (MATRIX_DENSE == right.m_type) {
            for (size_t j = 0; j < right.m_ncols; j++) {
                result.set(i, j, get(i, i) * right.get(i, j));
            }
        } else if (MATRIX_DIAGONAL == right.m_type) {
            result.set(i, i, m_data[i] * right.m_data[i]);
        } else if (MATRIX_LOWERTR == right.m_type) {
            for (size_t j = 0; j <= i; j++) {

                result.set(i, j, get(i, i) * right.get(i, j));
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
        Matrix result;
        if (isColumnVector() && right.isColumnVector()) { /* Sparse-sparse dot product */
            result = Matrix(1, 1, Matrix::MATRIX_SPARSE);
        } else {
            result = Matrix(m_nrows, right.m_ncols, Matrix::MATRIX_SPARSE);
        }
        result.m_sparse = r;
        result._createTriplet();
        return result;
    } else {
        // RHS is dense
        Matrix result(m_nrows, right.m_ncols);

        if (m_triplet != NULL && m_sparse == NULL)
            _createSparse();

        double alpha[2] = {1.0, 0.0};
        double beta[2] = {0.0, 0.0};

        if (right.m_dense == NULL) {
            right.m_dense = cholmod_allocate_dense(right.m_nrows, right.m_ncols, right.m_nrows, CHOLMOD_REAL, Matrix::cholmod_handle());
            for (size_t k = 0; k < right.length(); k++) {
                ((double*) right.m_dense->x)[k] = right.m_data[k];
            }
        }
        bool dotProd = isColumnVector() && right.isColumnVector();
        result.m_dense = cholmod_allocate_dense(
                dotProd ? 1 : result.m_nrows,
                dotProd ? 1 : result.m_ncols,
                dotProd ? 1 : result.m_nrows,
                CHOLMOD_REAL,
                Matrix::cholmod_handle());
        cholmod_sdmult(
                m_sparse,
                dotProd ? true : m_transpose,
                alpha,
                beta,
                right.m_dense,
                result.m_dense,
                Matrix::cholmod_handle()
                );
        for (size_t k = 0; k < result.length(); k++) {
            result.m_data[k] = (static_cast<double*> (result.m_dense->x))[k];
        }
        return result;
    }
}

bool Matrix::indexWithinBounds(size_t i, size_t j) {

    return (i >= 0 && j >= 0 && i < m_nrows && j < m_ncols) && !(m_type == MATRIX_LOWERTR && i < j);
}

Matrix::MatrixType Matrix::getType() const {

    return m_type;
}

void Matrix::init(size_t nr, size_t nc, MatrixType mType) {
    this -> m_transpose = false;
    this -> m_ncols = nc;
    this -> m_nrows = nr;
    this -> m_type = mType;
    switch (m_type) {
        case MATRIX_DENSE:
            m_dataLength = nc * nr;
            m_data = new double[m_dataLength]();
            break;
        case MATRIX_DIAGONAL:
            if (nc != nr) {
                throw std::invalid_argument("Diagonal matrices must be square!!!");
            }
            m_dataLength = nc;
            m_data = new double[m_dataLength]();
            break;
        case MATRIX_LOWERTR:
        case MATRIX_SYMMETRIC:
            if (nc != nr) {
                throw std::invalid_argument("Lower triangular and symmetric matrices must be square!!!");
            }
            m_dataLength = nc * (nc + 1) / 2;
            m_data = new double[m_dataLength]();

            break;
        case MATRIX_SPARSE:
            break;
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

inline void Matrix::_createTriplet() {
    _createSparse();
    if (m_sparse != NULL) { /* make triplets from sparse */
        m_triplet = cholmod_sparse_to_triplet(m_sparse, Matrix::cholmod_handle());
    }
}

bool Matrix::isSymmetric() {
    return (Matrix::MATRIX_SYMMETRIC == m_type)
            || (m_triplet != NULL && m_triplet->stype != 0)
            || (Matrix::MATRIX_DIAGONAL == m_type);
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
            ((double*) obj.m_triplet->x)[k] *= alpha;
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

}

Matrix Matrix::submatrixCopy(size_t row_start, size_t row_end, size_t col_start, size_t col_end) {
    /* DONE: Fully tested */
    if (row_end < row_start || col_end < col_start) {
        throw std::invalid_argument("Matrix::submatrixCopy:: start > end is not allowed");
    }
    if (row_end > m_nrows) {
        throw std::invalid_argument("Matrix::submatrixCopy:: row_end > total number of rows");
    }
    if (col_end > m_ncols) {
        throw std::invalid_argument("Matrix::submatrixCopy:: col_end > total number of columns");
    }
    
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
        dlacpy_((char*) "A",
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
        for (int i = 0; i <= rows; i++) {
            rs[i] = row_start + i;
        }
        for (int j = 0; j <= cols; j++) {
            cs[j] = col_start + j;
        }
        this->_createSparse();
        int nnz = std::max(1, (int) (rows * cols / 20));
        cholmod_sparse *sp = cholmod_allocate_sparse(rows, cols, nnz, 1, 1, m_sparse->stype, m_sparse->xtype, cholmod_handle());
        sp = cholmod_submatrix(m_sparse, rs, rows, cs, cols, 1, 1, Matrix::cholmod_handle());
        M.m_sparse = sp;
        M._createTriplet();
    } else {
        throw std::logic_error("Matrix::submatrixCopy is available only for MATRIX_DENSE and MATRIX_SPARSE type matrices.");
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

        //TODO tests!!!
        //TODO transposes!!!

        size_t left_cols = left_col_end - left_col_start + 1;
        size_t left_rows = left_row_end - left_row_start + 1;

        size_t right_cols = right_col_end - right_col_start + 1;
        size_t right_rows = right_row_end - right_row_start + 1;

        if (left_cols != right_rows) {
            throw std::invalid_argument("Dimension of sub-matrix mismatch (left_cols!=right_rows)");
        }

        size_t left_start_idx = left_row_start + left_col_start * m_nrows;
        size_t right_start_idx = right_row_start + right_col_start * right.m_nrows;

        Matrix result(left_rows, right_cols, MATRIX_DENSE);



        /*
         * C := beta*C + alpha*A*B
         * 
         * void cblas_dgemm(
         *  const enum CBLAS_ORDER Order,           CblasColMajor as always here
         *  const enum CBLAS_TRANSPOSE TransA,      is A traspose?
         *  const enum CBLAS_TRANSPOSE TransB,      is B traspose?
         *  const int M,                            # rows of matrix C
         *  const int N,                            # columns of matrix C
         *  const int K,                            # number of rows of matrix op(B)
         *  const double alpha,                     scalar alpha
         *  const double *A,                        matrix A
         *  const int lda,                          leading dimension of A
         *  const double *B,                        matrix B
         *  const int ldb,                          leading dimension of B
         *  const double beta,                      scalar beta (multiplies C)
         *  double *C,                              output matrix C
         *  const int ldc                           leading dimension of C
         * );
         */

        cblas_dgemm(
                CblasColMajor, CblasNoTrans, CblasNoTrans,
                left_rows,
                right_cols,
                left_cols,
                1.0f,
                m_data + left_start_idx,
                m_nrows,
                right.m_data + right_start_idx,
                right.m_nrows,
                0.0f,
                result.m_data,
                left_rows);
        return result;

    } else {
        throw std::logic_error("Matrix::multiplySubmatrix is implemented only for dense matrix multiplication.");
    }
}
