/* 
 * File:   Matrix.cpp
 * Author: Chung
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

/********* CONSTRUCTORS ************/
Matrix::Matrix() {
    m_nrows = 0;
    m_ncols = 0;
    m_data = new float;
    *m_data = 0;
    m_type = MATRIX_DENSE;
    m_dataLength = 0;
    m_transpose = false;
}

Matrix::Matrix(int nr, int nc) {
    init(nr, nc, MATRIX_DENSE);
}

Matrix::Matrix(int nr, int nc, MatrixType mType) {
    init(nr, nc, mType);
}

Matrix::Matrix(int nr, int nc, const float * dat) {
    init(nr, nc, MATRIX_DENSE);
    for (int j = 0; j < nc * nr; j++) {
        m_data[j] = dat[j];
    }
}

Matrix::Matrix(int nr, int nc, const float * dat, MatrixType mType) {
    init(nr, nc, mType);
    for (int j = 0; j < length(); j++) {
        m_data[j] = dat[j];
    }
}

Matrix::Matrix(const Matrix& orig) {
    m_ncols = orig.m_ncols;
    m_nrows = orig.m_nrows;
    m_transpose = orig.m_transpose;
    if (orig.m_type != MATRIX_SPARSE) {
        int n = orig.m_dataLength;
        if (n <= 0) {
            n = 1;
        }
        m_data = new float[n];
        for (int i = 0; i < n; i++) {
            m_data[i] = orig.m_data[i];
        }
        m_dataLength = orig.m_dataLength;
    }
    m_type = orig.m_type;
}

/********* DENSTRUCTOR ************/
Matrix::~Matrix() {
    this -> m_ncols = -1;
    this -> m_nrows = -1;
    if (m_data != NULL) {
        delete[] m_data;
    }
    m_data = NULL;
}

/********* GETTERS/SETTERS ************/
int Matrix::getNcols() const {
    return m_ncols;
}

int Matrix::getNrows() const {
    return m_nrows;
}

float * const Matrix::getData() const {
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

int Matrix::length() const {
    return m_dataLength;
}

/********* OTHER METHODS ************/

void Matrix::transpose() {
    if (m_type == MATRIX_DIAGONAL || m_type == MATRIX_SYMMETRIC) {
        return;
    }
    if (this -> m_transpose) {
        this -> m_transpose = false;
    } else {
        this -> m_transpose = true;
    }
    std::swap(this -> m_ncols, this -> m_nrows);
}

int Matrix::reshape(int nrows, int ncols) {
    int new_size = nrows*ncols;
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

float Matrix::get(const int i, const int j) const {
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
        int i_ = std::max(i, j);
        int j_ = std::min(i, j);
        return m_data[i_ + m_nrows * j_ - j_ * (j_ + 1) / 2];
    } else if (m_type == MATRIX_LOWERTR) {
        return m_transpose
                ? (j >= i) ? m_data[j + m_ncols * i - i * (i + 1) / 2] : 0.0f
                : (i >= j) ? m_data[i + m_nrows * j - j * (j + 1) / 2] : 0.0f;
    } else if (m_type == MATRIX_SPARSE) {
        if (m_sparseStorageType == CHOLMOD_TYPE_TRIPLET) {
            float val = 0.0f;
            for (int k = 0; k < m_triplet->nnz; k++) {
                if (i == ((int*) m_triplet->i)[k]) {
                    if (j == ((int*) m_triplet->j)[k]) {
                        val = ((double*) m_triplet->x)[k];
                        break;
                    }
                }
            }
            return val;
        } else {
            throw std::logic_error("Sorry... not supported yet!");
        }

    }
}

void Matrix::set(int i, int j, float v) {
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
        if (m_triplet->nnz == m_triplet->nzmax) {
            std::ostringstream oss;
            oss << "Cannot add more elements to this matrix; nnz_max = " << m_triplet->nzmax;
            throw std::out_of_range(oss.str());
        }
        ((int*) m_triplet->i)[m_triplet->nnz] = i;
        ((int*) m_triplet->j)[m_triplet->nnz] = j;
        ((double*) m_triplet->x)[m_triplet->nnz] = v;
        (m_triplet->nnz)++;
    } else {
        throw std::invalid_argument("Illegal operation");
    }

}

float Matrix::quad(Matrix& x) {
    if (!x.isColumnVector()) {
        throw std::invalid_argument("Method `quadratic` can only be applied to vectors!");
    }
    if (getNcols() != getNrows()) {
        throw std::invalid_argument("Method `quadratic` can only be applied on square matrices!");
    }
    if (x.getNrows() != m_ncols) {
        throw std::invalid_argument("The argument of quad(Matrix&) is not of appropriate dimension.");
    }
    float result = 0.0;
    if (MATRIX_DENSE == m_type || MATRIX_LOWERTR == m_type) {
        for (int i = 0; i < m_nrows; i++) {
            for (int j = 0; j < m_ncols; j++) {
                result += x[i] * get(i, j) * x[j];
            }
        }
    } else if (MATRIX_DIAGONAL == m_type) {
        for (int i = 0; i < m_nrows; i++) {
            result += x[i] * x[i] * get(i, i);
        }
    } else if (MATRIX_SYMMETRIC == m_type) {
        for (int i = 0; i < m_nrows; i++) {
            for (int j = 0; j < i; j++) {
                result += 2.0f * x[i] * get(i, j) * x[j];
            }
            result += x[i] * get(i, i) * x[i];
        }
    }
    return result;
}

float Matrix::quad(Matrix& x, Matrix& q) {
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
    float t = 0.0f;
    Matrix r;
    r = q*x;
    t = r[0] + quad(x);
    return t;
}

int Matrix::cholesky(Matrix& L) {
    if (m_type == MATRIX_SPARSE) {
        // Cholesky decomposition of a sparse matrix:
        m_sparse = cholmod_triplet_to_sparse(m_triplet, m_triplet->nzmax, m_cholmod_common);
        m_cholesky_factor = cholmod_analyze(m_sparse, m_cholmod_common); // analyze
        int status = cholmod_factorize(m_sparse, m_cholesky_factor, m_cholmod_common); // factorize
        L = *this;
        L.m_sparseStorageType = CHOLMOD_TYPE_FACTOR;
        return status;
    } else {
        L = *this;
        if (m_nrows != m_ncols) {
            throw std::invalid_argument("Method `cholesky` can only be applied to square matrices!");
        }
        int info;
#ifdef USE_LIBS
        if (m_type == MATRIX_DENSE) {
#ifdef DEBUG_VERBOSE
            std::cerr << "[Warning] Applying Cholesky on a possibly non-symmetric matrix!" << std::endl;
#endif
            info = LAPACKE_spotrf(LAPACK_COL_MAJOR, 'L', m_nrows, L.m_data, m_nrows);
            for (int i = 0; i < m_nrows; i++) {
                for (int j = 0; j > i; j++) {
                    L.set(i, j, 0.0f);
                }
            }
        } else if (m_type == MATRIX_SYMMETRIC) {
            info = LAPACKE_spptrf(LAPACK_COL_MAJOR, 'L', m_nrows, L.m_data);
            L.m_type = MATRIX_LOWERTR;
        }
#else
        //TODO: Make fail-safe (return proper status code!)
        info = 0;
        for (int i = 0; i < m_nrows; i++) {
            for (int j = 0; j < (i + 1); j++) {
                double s = 0;
                for (int k = 0; k < j; k++)
                    s += L[k * m_nrows + i] * L[k * m_nrows + j];
                L[j * m_nrows + i] = (i == j) ?
                        std::sqrt(m_data[i * m_nrows + i] - s) :
                        (1.0 / L[j * m_nrows + j] * (m_data[i * m_nrows + j] - s));
            }
        }
#endif
        return info;
    }
}

int Matrix::solveCholeskySystem(Matrix& solution, const Matrix & rhs) const {
    if (getType() == MATRIX_SPARSE) {
        cholmod_dense *x;
        cholmod_dense *b;
        b = cholmod_allocate_dense(rhs.m_nrows, rhs.m_ncols, rhs.m_nrows, CHOLMOD_REAL, m_cholmod_common);
        for (int k = 0; k < rhs.length(); k++) {
            ((double*) b->x)[k] = rhs[k];
        }
        x = cholmod_solve(CHOLMOD_A, m_cholesky_factor, b, m_cholmod_common);
        solution = Matrix(rhs.m_nrows, rhs.m_ncols);
        for (int k = 0; k < x->nzmax; k++) {
            solution.m_data[k] = ((double*) x->x)[k];
        }
        cholmod_free_dense(&x, m_cholmod_common);
        return 0;
    } else {
        int info = 1;
        solution = Matrix(rhs);
        if (m_nrows != m_ncols) {
            throw std::invalid_argument("Method `solveCholesky` can only be applied to square matrices!");
        }
#ifdef USE_LIBS
        info = LAPACKE_spotrs(LAPACK_COL_MAJOR, 'L', m_nrows, rhs.m_ncols, m_data, m_nrows, solution.m_data, m_nrows);
#else
        // Custom implementation!
#endif  
        return info;
    }
}

/********* OPERATORS ************/
bool Matrix::operator==(const Matrix& right) const {
    bool result = (m_type == right.m_type) &&
            (m_ncols == right.m_ncols) &&
            (m_nrows == right.m_nrows);
    for (unsigned int i = 0; i < m_nrows; i++) {
        for (unsigned int j = 0; j < m_ncols; j++) {
            result = result && (get(i, j) == right.get(i, j));
        }
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, const Matrix & obj) {
    os << "Matrix " << obj.m_nrows << "x" << obj.m_ncols << std::endl;
    if (obj.m_transpose) {
        os << "Stored as transpose : YES\n";
    }
    const char * const types[] = {"Dense", "Sparse", "Diagonal", "Lower Triangular", "Symmetric"};
    os << types[obj.m_type] << std::endl;
    for (int i = 0; i < obj.m_nrows; i++) {
        for (int j = 0; j < obj.m_ncols; j++) {
            os << std::setw(8) << std::setprecision(4) << obj.get(i, j) << ",";
        }
        os << std::endl;
    }
    return os;
}

float &Matrix::operator[](int sub) const {
    if (sub < 0 || sub > length()) {
        throw std::out_of_range("Exception: Index out of range for Matrix");
    }
    return m_data[sub];
}

Matrix & Matrix::operator+=(const Matrix & right) {
    if (m_ncols != right.m_ncols || m_nrows != right.m_nrows) {
        throw std::invalid_argument("Incompatible dimensions while using +=!");
    }
#ifdef USE_LIBS    
    cblas_saxpy(length(), 1.0f, right.m_data, 1, m_data, 1); // data = data + right.data
#else
    for (int i = 0; i < length(); i++) {
        m_data[i] += right[i];
    }
#endif
    return *this;
}

Matrix & Matrix::operator-=(const Matrix & right) {
    if (m_ncols != right.m_ncols || m_nrows != right.m_nrows) {
        throw std::invalid_argument("Incompatible dimensions while using +=!");
    }
#ifdef USE_LIBS    
    cblas_saxpy(length(), -1.0f, right.m_data, 1, m_data, 1); // data = data + right.data
#else
    for (int i = 0; i < length(); i++) {
        m_data[i] -= right[i];
    }
#endif
    return *this;
}

Matrix Matrix::operator+(const Matrix & right) const {
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
    float t = 0.0f;
    if (isColumnVector() && right.isColumnVector() && length() == right.length()) {
        // multiplication of two column vectors = dot product
        Matrix r(1, 1);
        for (int i = 0; i < m_nrows * m_ncols; i++) {
            t += m_data[i] * right.m_data[i];
        }
        r[0] = t;
        return r;
    }
    if (m_ncols != right.m_nrows) {
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
    if (m_type != MATRIX_SPARSE) {
        m_dataLength = right.m_dataLength;
        if (m_data != NULL) {
            delete m_data;
        }
        m_data = new float[m_dataLength];
    }
    m_transpose = right.m_transpose;
    m_triplet = right.m_triplet;
    m_sparse = right.m_sparse;
    m_sparseStorageType = right.m_sparseStorageType;
    m_cholmod_common = right.m_cholmod_common;
    m_cholesky_factor = right.m_cholesky_factor;

#ifdef USE_LIBS
    cblas_scopy(m_dataLength, right.m_data, 1, m_data, 1);
#else
    for (int i = 0; i < m_dataLength; i++) {
        m_data[i] = right[i];
    }
#endif
    return *this;
}

/********* PRIVATE METHODS ************/
void Matrix::domm(const Matrix &right, Matrix &result) const {
    // multiply with LHS being dense
    float t;
    //std::cout << "multiplication - left = " << m_nrows << "x" << m_ncols << std::endl;
    for (int i = 0; i < m_nrows; i++) {
        for (int j = 0; j < right.m_ncols; j++) {
            t = 0.0;
            //std::cout << "For element C(" << i << "," << j << ")\n";
            for (int k = 0; k < m_ncols; k++) {
                //std::cout << "  += A(" << i << "," << k << ") * B(" << k << "," << i << ")" << std::endl;
                if (!(right.getType() == MATRIX_LOWERTR && k < j)) {
                    t += get(i, k) * right.get(k, j);
                }
            }
            result.set(i, j, t);
        }
    }
}

Matrix Matrix::multiplyLeftDense(const Matrix & right) const{
    if (MATRIX_DENSE == right.m_type) { // RHS is also dense
        Matrix result(m_nrows, right.m_ncols);
#ifdef USE_LIBS
        cblas_sgemm(CblasColMajor,
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
        for (int j = 0; j < m_ncols; j++) {
            for (int i = 0; i < m_nrows; i++) {
                result.set(i, j, result.get(i, j) * right.get(j, j));
            }
        }
        return result;
    } else if (MATRIX_SYMMETRIC == right.m_type || MATRIX_LOWERTR == right.m_type) {
        Matrix result(m_nrows, right.m_ncols, Matrix::MATRIX_DENSE);
        domm(right, result);
        return result;
    } else {
        throw std::logic_error("Sorry... this is not implemented yet!");
    }
}

Matrix Matrix::multiplyLeftSymmetric(const Matrix& right) const {
    // multiply when the LHS is symmetric    
    Matrix result(m_nrows, right.m_ncols);
    if (right.isColumnVector()) {
#ifdef USE_LIBS
        cblas_sspmv(CblasColMajor, CblasLower,
                m_nrows, 1.0f, m_data,
                right.m_data, 1,
                0.0, result.m_data, 1);
        return result;
#endif
    }
    domm(right, result);
    return result;
}

Matrix Matrix::multiplyLeftDiagonal(const Matrix & right) const{
    // multiply when the LHS is diagonal
    Matrix result(m_nrows, right.m_ncols, right.m_type);
    for (int i = 0; i < m_nrows; i++) {
        if (MATRIX_SYMMETRIC == right.m_type) {
            for (int j = i; j < right.m_ncols; j++) {
                result.set(i, j, get(i, i) * right.get(i, j));
            }
        } else if (MATRIX_DENSE == right.m_type) {
            for (int j = 0; j < right.m_ncols; j++) {
                result.set(i, j, get(i, i) * right.get(i, j));
            }
        } else if (MATRIX_DIAGONAL == right.m_type) {
            result.set(i, i, m_data[i] * right.m_data[i]);
        } else if (MATRIX_LOWERTR == right.m_type) {
            for (int j = 0; j <= i; j++) {
                result.set(i, j, get(i, i) * right.get(i, j));
            }
        }
    }
    return result;
}

Matrix Matrix::multiplyLeftSparse(Matrix& right) {
    if (right.m_type == MATRIX_SPARSE) {
        // RHS is sparse
        throw std::logic_error("S*S not implemented yet!");
    } else {
        // RHS is dense
        Matrix result(m_nrows, right.m_ncols);

        if (m_triplet != NULL && m_sparse == NULL)
            createSparseFromTriplet();

        double alpha[2] = {1.0, 0.0};
        double beta[2] = {0.0, 0.0};

        if (right.m_dense == NULL) {
            right.m_dense = cholmod_allocate_dense(right.m_nrows, right.m_ncols, right.m_nrows, CHOLMOD_REAL, m_cholmod_common);
            for (int k = 0; k < right.length(); k++) {
                ((double*) right.m_dense->x)[k] = right.m_data[k];
            }
        }
        result.m_dense = cholmod_allocate_dense(result.m_nrows, result.m_ncols, result.m_nrows, CHOLMOD_REAL, m_cholmod_common);
        cholmod_sdmult(
                m_sparse,
                m_transpose ? 1 : 0,
                alpha,
                beta,
                right.m_dense,
                result.m_dense,
                m_cholmod_common
                );
        for (int k = 0; k < result.length(); k++) {
            result.m_data[k] = ((double*)result.m_dense->x)[k];
        }
        return result;
    }
}

bool Matrix::indexWithinBounds(int i, int j) {
    return (i >= 0 && j >= 0 && i < m_nrows && j < m_ncols) && !(m_type == MATRIX_LOWERTR && i < j);
}

Matrix::MatrixType Matrix::getType() const {
    return m_type;
}

void Matrix::init(int nr, int nc, MatrixType mType) {
    this -> m_transpose = false;
    this -> m_ncols = nc;
    this -> m_nrows = nr;
    this -> m_type = mType;
    switch (m_type) {
        case MATRIX_DENSE:
            m_dataLength = nc * nr;
            m_data = new float[m_dataLength]();
            break;
        case MATRIX_DIAGONAL:
            if (nc != nr) {
                throw std::invalid_argument("Diagonal matrices must be square!!!");
            }
            m_dataLength = nc;
            m_data = new float[m_dataLength]();
            break;
        case MATRIX_LOWERTR:
        case MATRIX_SYMMETRIC:
            if (nc != nr) {
                throw std::invalid_argument("Lower triangular matrices must be square!!!");
            }
            m_dataLength = nc * (nc + 1) / 2;
            m_data = new float[m_dataLength]();
            break;
        case MATRIX_SPARSE:
            break;
    }

}

void Matrix::createSparseFromTriplet() {
    if (m_type == MATRIX_SPARSE && m_triplet != NULL) {
        m_sparse = cholmod_triplet_to_sparse(m_triplet, m_triplet->nzmax, m_cholmod_common);
    }
}
