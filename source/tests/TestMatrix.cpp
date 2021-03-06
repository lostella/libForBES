/*
 * File:   TestMatrix.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Sep 14, 2015, 3:34:41 PM
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

#include "TestMatrix.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestMatrix);

void _testSubmatrixMultiply(Matrix& A, Matrix& B);

TestMatrix::TestMatrix() {
}

TestMatrix::~TestMatrix() {
}

void TestMatrix::setUp() {
}

void TestMatrix::tearDown() {
    Matrix::destroy_handle();
}

void TestMatrix::testToggleDiagonal() {
    // Diagonal to vector
    size_t n = 10;
    double alpha = 3.0;
    Matrix X = MatrixFactory::MakeIdentity(n, alpha);
    _ASSERT_EQ(Matrix::MATRIX_DIAGONAL, X.getType());
    X.toggle_diagonal();
    _ASSERT_EQ(Matrix::MATRIX_DENSE, X.getType());
    _ASSERT(X.isColumnVector());
    _ASSERT_EQ(n, X.getNrows());
    for (size_t i = 0; i < n; i++) {
        _ASSERT_EQ(alpha, X[i]);
    }

    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, x.getType());
    _ASSERT(x.isColumnVector());
    x.toggle_diagonal();
    _ASSERT_EQ(Matrix::MATRIX_DIAGONAL, x.getType());
    x.toggle_diagonal();
    _ASSERT(x.isColumnVector());


    Matrix y = MatrixFactory::MakeRandomMatrix(2, n, 0.0, 1.0);
    _ASSERT_EXCEPTION(y.toggle_diagonal(), std::invalid_argument);

    Matrix z = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_LOWERTR);
    _ASSERT_EXCEPTION(z.toggle_diagonal(), std::invalid_argument);


}

void TestMatrix::testGetSet() {
    const size_t n = 10;
    Matrix f(n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < 10; j++) {
            f.set(i, j, static_cast<double> (3 * i + 5 * j + 13));
        }
    }
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(static_cast<double> (3 * i + 5 * j + 13), f.get(i, j));
        }
    }

    Matrix o;
    _ASSERT(o.isEmpty());
    _ASSERT_EXCEPTION(o.get(0, 0), std::out_of_range);
}

void TestMatrix::testOpplus() {
    size_t n = 1000;
    size_t m = 50;
    Matrix x = MatrixFactory::MakeRandomMatrix(n, m, -1.5, 2.0);
    Matrix y(x);

    x.plusop();

    const double tol = 1e-12;
    for (size_t i = 0; i < x.getNrows(); i++) {
        for (size_t j = 0; j < x.getNcols(); j++) {
            if (y.get(i, j) >= 0.0) {
                _ASSERT_NUM_EQ(y.get(i, j), x.get(i, j), tol);
            } else {
                _ASSERT_NUM_EQ(0.0, x.get(i, j), tol);
            }
        }
    }
}

void TestMatrix::testOpplus2() {
    size_t n = 10;
    size_t m = 1000;
    Matrix x = MatrixFactory::MakeRandomMatrix(n, m, -1.5, 2.0);   
    
    Matrix z(n, m);
    x.plusop(&z);
    
    const double tol = 1e-12;
    for (size_t i = 0; i < x.getNrows(); i++) {
        for (size_t j = 0; j < x.getNcols(); j++) {
            if (x.get(i, j) >= 0.0) {
                _ASSERT_NUM_EQ(x.get(i, j), z.get(i, j), tol);
            } else {
                _ASSERT_NUM_EQ(0.0, z.get(i, j), tol);
            }
        }
    }
}

void TestMatrix::testOpplusSparse() {
    size_t n = 1000;
    size_t m = 50;
    size_t nnz = 200;
    Matrix x = MatrixFactory::MakeRandomSparse(n, m, nnz, -1.5, 2.0);
    Matrix y(x);

    x.plusop();

    const double tol = 1e-12;
    for (size_t i = 0; i < x.getNrows(); i++) {
        for (size_t j = 0; j < x.getNcols(); j++) {
            if (y.get(i, j) >= 0.0) {
                _ASSERT_NUM_EQ(y.get(i, j), x.get(i, j), tol);
            } else {
                _ASSERT_NUM_EQ(0.0, x.get(i, j), tol);
            }
        }
    }
}

void TestMatrix::testQuadratic() {
    /* Test quadratic with diagonal matrices */
    Matrix *f;    
    for (size_t n = 5; n < 12; n++) {
	Matrix *x = new Matrix(n, 1);
        f = new Matrix(n, n);
        for (size_t i = 0; i < n; i++) {
            _ASSERT_OK(f -> set(i, i, 1.0));
            (*x)[i] = i + 1;
        }
        double r = f -> quad(*x);
        _ASSERT_EQ(static_cast<double> (n * (n + 1)* (2 * n + 1) / 6), 2 * r);
        _ASSERT_OK(delete f);
        _ASSERT_OK(delete x);
    }
    Matrix A(5, 6);
    Matrix y(6, 1);
    _ASSERT_EXCEPTION(A.quad(y), std::invalid_argument);

    Matrix B(5, 5);
    Matrix z(5, 2);
    _ASSERT_EXCEPTION(B.quad(z), std::invalid_argument);
    _ASSERT_EXCEPTION(B.quad(z, z), std::invalid_argument);

    Matrix C(5, 5);
    Matrix s(7, 1);
    _ASSERT_EXCEPTION(C.quad(s), std::invalid_argument);

}

void TestMatrix::testQuadratic2() {
    double fdata[9] = {-1.0, 3.0, 1.0, 2.0, -1.0, -1.0, 5.0, 2.0, -5.0};
    double xdata[3] = {1.0, 2.0, 3.0};

    Matrix f(3, 3, fdata);
    Matrix x(3, 1, xdata);

    double r;
    _ASSERT_OK(r = f.quad(x));
    _ASSERT_EQ(-8.0, r);
}

void TestMatrix::testQuadratic3() {
    double fdata[9] = {-1.0, -3.0, 7.5, 2.0, -1.0, -1.0, 5.0, 2.0, -5.0};
    double xdata[3] = {-1.5, 2.0, 3.0};
    double qdata[3] = {5.0, -6.0, 13.5};

    Matrix f(3, 3, fdata);
    Matrix x(3, 1, xdata);
    Matrix q(3, 1, qdata);


    double r;
    _ASSERT_OK(r = f.quad(x, q));
    _ASSERT_EQ(-28.25, r);
}

void TestMatrix::testQuadraticDot() {
    double ydata[4] = {-2.0, 5.0, -6.0, -11.0};
    double xdata[4] = {10.0, 2.0, 3.0, 4.0};

    const size_t n = 4;
    const size_t m = 1;

    Matrix y(n, m, ydata);
    Matrix x(n, m, xdata);

    Matrix r = y*x;

    _ASSERT_EQ(m, r.getNcols());
    _ASSERT_EQ(m, r.getNrows());
    _ASSERT_EQ(m, r.length());
    _ASSERT_EQ(-72.0, r[0]);
}

void TestMatrix::test_MDD1() {
    double fdata[9] = {-1.0, 3.0, 1.0, 2.0, -1.0, -1.0, 5.0, 2.0, -5.0};
    double xdata[3] = {1.0, 2.0, 3.0};

    const size_t n = 3;
    const size_t m = 1;

    Matrix f(n, n, fdata);
    Matrix x(n, m, xdata);

    Matrix r;
    r = f*x;

    _ASSERT_EQ(n, r.getNrows());
    _ASSERT_EQ(m, r.getNcols());
    _ASSERT_EQ(n, r.length());
    _ASSERT_EQ(18.0, r[0]);
    _ASSERT_EQ(7.0, r[1]);
    _ASSERT_EQ(-16.0, r[2]);

    Matrix A(5, 6);
    Matrix B(3, 5);
    Matrix C;
    _ASSERT_EXCEPTION(C = A*B, std::invalid_argument);
}

void TestMatrix::testAssignment() {
    const size_t nRows = 3;
    const size_t nCols = 2;
    const size_t size = nRows * nCols;

    Matrix f = MatrixFactory::MakeRandomMatrix(nRows, nCols, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix g;
    g = f;

    _ASSERT_EQ(Matrix::MATRIX_DENSE, g.getType());
    _ASSERT_EQ(size, g.length());
    for (size_t i = 0; i < size; i++) {
        _ASSERT(g[i] >= 0);
        _ASSERT(g[i] <= 1);
    }
    _ASSERT_OK(f = f);
}

void TestMatrix::testAdditionBad() {
    Matrix A(5, 6);
    Matrix B(7, 8);
    Matrix C;
    CPPUNIT_ASSERT_THROW(C = A + B, std::invalid_argument);
    CPPUNIT_ASSERT_THROW(C = A - B, std::invalid_argument);
    CPPUNIT_ASSERT_THROW(A += B, std::invalid_argument);
    CPPUNIT_ASSERT_THROW(A -= B, std::invalid_argument);
}

void TestMatrix::test_ADD1() {
    const size_t nRows = 3;
    const size_t nCols = 2;
    const size_t size = nRows * nCols;
    double *a, *b;
    a = new double[size];
    b = new double[size];
    for (size_t i = 0; i < size; i++) {
        a[i] = i;
        b[i] = 3 * i + 7;
    }

    Matrix A(nRows, nCols, a);
    Matrix B(nRows, nCols, b);


    Matrix C;
    C = A + B;

    for (size_t i = 0; i < size; i++) {
        _ASSERT_EQ(a[i] + b[i], C[i]);
    }

    Matrix D = A + B + C;
    for (size_t i = 0; i < size; i++) {
        _ASSERT_EQ(2 * C[i], D[i]);
    }
    delete[] a;
    delete[] b;
}

void TestMatrix::testFBMatrix() {
    /* Test FBMatrix() - empty constructor */
    Matrix *fBMatrix = new Matrix();
    _ASSERT_EQ(static_cast<size_t> (0), fBMatrix->getNcols());
    _ASSERT_EQ(static_cast<size_t> (0), fBMatrix->getNrows());
    delete fBMatrix;

    /* Test FBMatrix(size_t, size_t, double*) - Provide data */
    const size_t n = 5;
    double *x = new double[n];
    for (size_t i = 0; i < n; i++) {
        x[i] = 1 + 7 * i;
    }
    Matrix f(n, 1, x);
    delete[] x;
    for (size_t i = 0; i < n; i++) {
        _ASSERT_EQ(static_cast<double> (1 + 7 * i), f[i]);
    }

    double s = 0.0;

    _ASSERT_EXCEPTION(fBMatrix = new Matrix(3, 4, Matrix::MATRIX_DIAGONAL), std::invalid_argument);
    _ASSERT_EXCEPTION(fBMatrix = new Matrix(3, 4, Matrix::MATRIX_SYMMETRIC), std::invalid_argument);
    _ASSERT_EXCEPTION(fBMatrix = new Matrix(3, 4, Matrix::MATRIX_LOWERTR), std::invalid_argument);
//    _ASSERT_EXCEPTION(s = f[-1], std::out_of_range);
//    _ASSERT_EXCEPTION(s = f[n], std::out_of_range);
    _ASSERT_NUM_EQ(0.0, s, 1e-10);
    _ASSERT_OK(Matrix::destroy_handle());
    _ASSERT_EQ(0, Matrix::destroy_handle());

    Matrix E;
    Matrix T(E);
    _ASSERT(T.isEmpty());
}

void TestMatrix::testMakeRandomFBMatrix() {
    const size_t nRows = 10;
    const size_t nCols = 20;
    const double offset = 0.1;
    const double scale = 3.5;
    Matrix f = MatrixFactory::MakeRandomMatrix(nRows, nCols, offset, scale, Matrix::MATRIX_DENSE);

    _ASSERT_EQ(nCols, f.getNcols());
    _ASSERT_EQ(nRows, f.getNrows());
    for (size_t i = 0; i < nRows * nCols; i++) {
        if (i > 0) {
            _ASSERT_NEQ(f[i], f[i - 1]);
        }
        _ASSERT(f[i] >= offset);
        _ASSERT(f[i] <= offset + scale);
    }
}

void TestMatrix::testGetData() {
    const size_t n = 100;
    double *myData = new double[n];
    for (size_t j = 0; j < n; j++) {
        myData[j] = j;
    }
    Matrix * mat = new Matrix(n, 1, myData);

    double * retrievedData;
    retrievedData = mat->getData();

    for (size_t j = 0; j < n; j++) {
        _ASSERT_EQ(j, static_cast<size_t> (retrievedData[j]));
    }

    retrievedData[0] = 666;
    _ASSERT_EQ(666, static_cast<int> ((*mat)[0]));

    delete(mat);
}

void TestMatrix::testGetNcols() {
    size_t y = static_cast<size_t> (5 + 50 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX));
    Matrix mat(10, y);
    _ASSERT_EQ(y, mat.getNcols());
}

void TestMatrix::testGetNrows() {
    size_t r = static_cast<size_t> (5 + 50 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX));
    Matrix mat(r, 10);
    _ASSERT_EQ(r, mat.getNrows());

    const size_t n = 5;
    Matrix *gm = new Matrix(n, n);
    _ASSERT_NOT(gm->isEmpty());
    _ASSERT_EQ(n, gm->getNrows());
    _ASSERT_EQ(n, gm->getNcols());
    delete gm;
}

void TestMatrix::testIsColumnVector() {
    Matrix rowVec(2345, 1);
    _ASSERT(rowVec.isColumnVector());
}

void TestMatrix::testIsRowVector() {
    Matrix rowVec(1, 100);
    _ASSERT(rowVec.isRowVector());
}

void TestMatrix::testLength() {
    size_t nRep = 20;    

    for (size_t i = 0; i < nRep; i++) {
	size_t x = static_cast<size_t> (5 + 50 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX));
        size_t y = static_cast<size_t> (5 + 50 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX));
        Matrix * f = new Matrix(x, y);
        _ASSERT_EQ(x*y, f->length());
        delete f;
    }

}

void TestMatrix::testReshape() {
    size_t n = 45;
    size_t m = 56;
    Matrix f(n, m);
    size_t k = 5;
    int status = f.reshape(k, k);
    _ASSERT_EQ(0, status);
    _ASSERT_EQ(k, f.getNcols());
}

void TestMatrix::testReshapeBad() {
    size_t n = 45;
    size_t m = 56;
    Matrix f(n, m);
    int status = f.reshape(n, m + 1);
    _ASSERT_EQ(-2, status);
    _ASSERT_EQ(n, f.getNrows());
    _ASSERT_EQ(m, f.getNcols());

    status = f.reshape(0, 1);
    _ASSERT_EQ(-1, status);
}

void TestMatrix::testDiagonalGetSet() {
    const size_t n = 10;
    double *myData = new double[n];
    for (size_t j = 0; j < n; j++) {
        myData[j] = j + 1.0;
    }

    Matrix *A = new Matrix(n, n, myData, Matrix::MATRIX_DIAGONAL);
    delete[] myData;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(static_cast<double> (i + 1)*(i == j), A -> get(i, j));
        }
    }

    _ASSERT_EQ(Matrix::MATRIX_DIAGONAL, A -> getType());

    double t = -1.234;
    _ASSERT_OK(A -> set(3, 3, t));
    _ASSERT_EQ(t, A -> get(3, 3));
    _ASSERT_EXCEPTION(A -> set(3, 4, 1.0), std::invalid_argument);
    _ASSERT_EXCEPTION(A -> get(n, n - 1), std::out_of_range);
    _ASSERT_EXCEPTION(A -> set(n, n, 100.0), std::out_of_range);

    delete A;
}

void TestMatrix::testDiagonalMultiplication() {
    /* Diag * Dense = Dense */
    const size_t n = 10;
    const size_t m = 3;
    double *myData = new double[n];
    for (size_t j = 0; j < n; j++) {
        myData[j] = j + 1.0;
    }
    Matrix *A = new Matrix(n, n, myData, Matrix::MATRIX_DIAGONAL);
    Matrix B(n, m);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            B.set(i, j, i + 1.0);
        }
    }

    Matrix C;
    _ASSERT_OK(C = (*A) * B);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, C.getType());
    _ASSERT_EQ(n, C.getNrows());
    _ASSERT_EQ(m, C.getNcols());
    _ASSERT_NOT(C.isEmpty());
    _ASSERT_EQ(n*m, C.length());

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            _ASSERT_EQ(std::pow(i + 1.0, 2.0), C.get(i, j));
        }
    }

    _ASSERT_OK(delete A);
    delete[] myData;
}

void TestMatrix::testDiagonalMultiplication2() {
    /* Diag * Diag = Diag */
    const size_t n = 10;
    double *aData = new double[n];
    double *bData = new double[n];

    for (size_t i = 0; i < n; i++) {
        aData[i] = i + 1.0;
        bData[i] = n - i;
    }
    Matrix *A = new Matrix(n, n, aData, Matrix::MATRIX_DIAGONAL);
    Matrix *B = new Matrix(n, n, bData, Matrix::MATRIX_DIAGONAL);

    delete[] aData;
    delete[] bData;

    Matrix C = (*A) * (*B);
    _ASSERT_EQ(Matrix::MATRIX_DIAGONAL, A -> getType());
    _ASSERT_EQ(Matrix::MATRIX_DIAGONAL, B -> getType());
    _ASSERT_EQ(Matrix::MATRIX_DIAGONAL, C.getType());

    for (size_t i = 0; i < n; i++) {
        _ASSERT_EQ((i + 1.0)*(n - i), C.get(i, i));
    }
}

void TestMatrix::testDenseTimesDiagonal() {
    const size_t nRows = 7;
    const size_t nCols = 3;
    const size_t size = nRows * nCols;
    double *aData;
    double *bData;

    aData = new double[size];
    for (size_t i = 0; i < size; i++) {
        aData[i] = i;
    }

    bData = new double[nCols];
    for (size_t i = 0; i < nCols; i++) {
        bData[i] = 3.0 * (i + 1);
    }


    Matrix A(nRows, nCols, aData, Matrix::MATRIX_DENSE);
    Matrix B(nCols, nCols, bData, Matrix::MATRIX_DIAGONAL);

    Matrix C = A * B;

    for (size_t i = 0; i < C.getNrows(); i++) {
        for (size_t j = 0; j < C.getNcols(); j++) {
            _ASSERT_EQ(A.get(i, j) * B.get(j, j), C.get(i, j));
        }
    }

    delete[] aData;
    delete[] bData;
}

void TestMatrix::testQuadDiagonal() {
    const size_t n = 10;
    double *aData = new double[n];
    double *xData = new double[n];
    for (size_t i = 0; i < n; i++) {
        aData[i] = i + 1.0;
        xData[i] = 3.0 * (i + 1.0);
    }
    Matrix *A = new Matrix(n, n, aData, Matrix::MATRIX_DIAGONAL);
    Matrix *x = new Matrix(n, 1, xData);
    double val = A -> quad(*x, *x);

    const double correct = 17077.5;
    _ASSERT_EQ(correct, val);

    val = A -> quad(*x);
    const double correct2 = 13612.5;
    _ASSERT_EQ(correct2, val);

    delete A;
    delete x;
    _ASSERT_OK(delete[] aData);
    _ASSERT_OK(delete[] xData);
}

void TestMatrix::testSubtract() {
    const size_t n = 10;
    Matrix X = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix Y = X - X;
    for (size_t i = 0; i < n * n; i++) {
        _ASSERT_EQ(0.0, Y.get(i));
    }

    X = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    X -= X;
    for (size_t i = 0; i < n * n; i++) {
        _ASSERT_EQ(0.0, X.get(i));
    }

    X = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix B = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix Z = X - B;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(X.get(i, j) - B.get(i, j), Z.get(i, j));
        }
    }

}

void TestMatrix::testLowerTriangular_getSet() {
    for (size_t n = 1; n < 20; n++) {
        Matrix *A = new Matrix(n, n, Matrix::MATRIX_LOWERTR);

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j <= i; j++) {
                A -> set(i, j, 3.2 * i + 7.5 * j + 0.45);
            }
        }

        const double tol = 1e-6;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                if (i >= j) {
                    _ASSERT_NUM_EQ(3.2 * i + 7.5 * j + 0.45, A->get(i, j), tol);
                } else {
                    _ASSERT_EQ(0.0, A->get(i, j));
                }
            }
        }

        size_t exp_length = n * (n + 1) / 2;
        _ASSERT_EQ(exp_length, A->length());
        delete A;
    }
}

void TestMatrix::testLowerTriangularTraspose_getSet() {
    const size_t n = 10;
    Matrix A = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_LOWERTR);
    Matrix AT(A);
    AT.transpose();
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(A.get(i, j), AT.get(j, i));
        }
    }
    AT.transpose();
    _ASSERT_EQ(A, AT);
}

void TestMatrix::testSymmetric_getSet() {
    const size_t n = 4;
    Matrix *A = new Matrix(n, n, Matrix::MATRIX_SYMMETRIC);
    _ASSERT(A->isSymmetric());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            A -> set(i, j, 3.2 * i + 0.2 * j + 0.45);
        }
    }
    const double tol = 1e-6;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i >= j) {
                _ASSERT_NUM_EQ(3.2f * i + 0.2f * j + 0.45f, A->get(i, j), tol);
            } else {
                _ASSERT_NUM_EQ(A -> get(j, i), A->get(i, j), tol);
            }
        }
    }
    size_t exp_length = n * (n + 1) / 2;
    _ASSERT_EQ(exp_length, A->length());

    delete A;
}

void TestMatrix::testTranspose() {
    size_t n = 5;
    size_t m = 6;
    Matrix A(n, m);
    A.transpose();
    _ASSERT_EQ(n, A.getNcols());
    _ASSERT_EQ(m, A.getNrows());
    A.transpose();
    _ASSERT_EQ(n, A.getNrows());
    _ASSERT_EQ(m, A.getNcols());

    Matrix X(n, m);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            X.set(i, j, 7.5 * i - 2.8 * j - 1.0);
        }
    }
    X.transpose();
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_OK(X.get(i, j));
        }
    }

    Matrix Y(n, m);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            Y.set(i, j, 7.5 * i - 2.8 * j - 1.0);
        }
    }

    Matrix YT(m, n); // Construct the transpose of Y as YT
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            YT.set(i, j, Y.get(j, i));
        }
    }

    Y.transpose();
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(YT.get(i, j), Y.get(i, j));
        }
    }

}

void TestMatrix::test_MXH() {
    const size_t n = 10;
    Matrix D(n, n, Matrix::MATRIX_DIAGONAL);
    for (size_t i = 0; i < n; i++) {
        D[i] = i + 1.0;
    }

    Matrix S(n, n, Matrix::MATRIX_SYMMETRIC);
    _ASSERT(S.isSymmetric());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            S.set(i, j, -3.1 * i + 3.25 * j + 5.35);
        }
    }

    Matrix DS = D * S;
    _ASSERT_EQ(Matrix::MATRIX_SYMMETRIC, DS.getType());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            _ASSERT_EQ((i + 1.0)*(-3.1 * i + 3.25 * j + 5.35), DS.get(i, j));
        }
    }
}

void TestMatrix::test_MDX() { // dense * diagonal
    const size_t n = 10;
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 2.0, Matrix::MATRIX_DENSE);
    Matrix X = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 2.0, Matrix::MATRIX_DIAGONAL);
    Matrix result = D*X;
    Matrix XX(n, n, Matrix::MATRIX_DENSE);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            XX.set(i, j, X.get(i, j));
        }
    }
    Matrix correct = D*XX;
    _ASSERT_EQ(correct, result);

    result = X*D;
    correct = XX*D;
    _ASSERT_EQ(correct, result);

    const size_t m = 12;
    D = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 2.0, Matrix::MATRIX_DENSE);
    X = MatrixFactory::MakeRandomMatrix(m, m, 0.0, 2.0, Matrix::MATRIX_DIAGONAL);
    XX = Matrix(m, m, Matrix::MATRIX_DENSE);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            XX.set(i, j, X.get(i, j));
        }
    }

    correct = D*XX;
    result = D*X;
    _ASSERT_EQ(correct, result);

    _ASSERT_EQ(n, D.getNrows());
    _ASSERT_EQ(m, D.getNcols());
    D.transpose();
    _ASSERT_EQ(m, D.getNrows());
    _ASSERT_EQ(n, D.getNcols()); /* Now D is m-by-n */
    X = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 2.0, Matrix::MATRIX_DIAGONAL); /* X is n-by-n */
    result = D*X;
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < n; i++) {
            _ASSERT_NUM_EQ(D.get(i, j) * X.get(j, j), result.get(i, j), 1e-8);
        }
    }
}

void TestMatrix::test_MSX() { /* Sparse * Diagonal */
    const size_t n = 10;
    Matrix X(n, n, Matrix::MATRIX_DIAGONAL);
    for (size_t i = 0; i < n; i++) {
        X.set(i, i, i + 1.0);
    }
    Matrix S = MatrixFactory::MakeRandomSparse(n, n, static_cast<size_t> (2.5 * n), 0.0, 1.0);
    Matrix R = S*X;
    _ASSERT_EQ(Matrix::MATRIX_SPARSE, R.getType());
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < n; i++) {
            _ASSERT_NUM_EQ(S.get(i, j) * X.get(j, j), R.get(i, j), 1e-8);
        }
    }
}

void TestMatrix::test_MXL() {
    size_t n = 6;
    Matrix D(n, n, Matrix::MATRIX_DIAGONAL);
    for (size_t i = 0; i < n; i++) {
        D[i] = i + 1.0;
    }

    Matrix L(n, n, Matrix::MATRIX_LOWERTR);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            L.set(i, j, -0.1 * i + 0.45 * j + 1.01);
        }
    }

    Matrix DL = D * L;
    _ASSERT_EQ(Matrix::MATRIX_LOWERTR, DL.getType());
    _ASSERT_EQ(n, DL.getNcols());
    _ASSERT_EQ(n, DL.getNrows());
    _ASSERT_NOT(DL.isEmpty());

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            _ASSERT_EQ(0.0, DL.get(i, j) - D[i] * L.get(i, j));
        }
    }
}

void TestMatrix::test_MDH() {
    const size_t n = 4;
    double a[n * n] = {5, 11, -2, -1,
        6, 3, 7, -1,
        -21, 0, 13, -1,
        -18, -1, -17, 30};
    Matrix A(n, n, a, Matrix::MATRIX_DENSE);

    Matrix S(n, n, Matrix::MATRIX_SYMMETRIC);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            S.set(i, j, 3 * i + 2 * j + 4.0);
        }
    }

    Matrix AS = A*S;
    _ASSERT_NOT(AS.isEmpty());
    _ASSERT_NOT(AS.isColumnVector());
    _ASSERT_NOT(AS.isRowVector());
    _ASSERT_EQ(Matrix::MATRIX_DENSE, AS.getType());
    _ASSERT_EQ(n, AS.getNcols());
    _ASSERT_EQ(n, AS.getNrows());

    double asData[n * n] = {-382.0, 52.0, -50.0, 369.0,
        -433.0, 89.0, -50.0, 422.0,
        -478.0, 129.0, -43.0, 474.0,
        -544.0, 169.0, -23.0, 525.0};
    Matrix AS_correct(n, n, asData, Matrix::MATRIX_DENSE);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(AS_correct.get(i, j), AS.get(i, j));
        }
    }
}

void TestMatrix::test_MDL() {
    const size_t n = 4;
    double a[n * n] = {5, 11, -2, -1,
        6, 3, 7, -1,
        -21, 0, 13, -1,
        -18, -1, -17, 30};
    Matrix A(n, n, a, Matrix::MATRIX_DENSE);

    Matrix L(n, n, Matrix::MATRIX_LOWERTR);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            L.set(i, j, 3 * i + 2 * j + 4.0);
        }
    }

    Matrix AL = A*L;
    double alCorrect[n * n] = {-382.0, 52.0, -50.0, 369.0, -468.0, 12.0, -36.0, 429.0,
        -600.0, -17.0, -107.0, 496.0, -342.0, -19.0, -323.0, 570.0};
    Matrix AL_correct(n, n, alCorrect, Matrix::MATRIX_DENSE);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            _ASSERT_EQ(AL.get(i, j), AL_correct.get(i, j));
        }
    }
}

void TestMatrix::testQuadSymmetric() {
    const size_t n = 10;

    Matrix A(n, n, Matrix::MATRIX_SYMMETRIC);
    Matrix x(n, 1);

    for (size_t i = 0; i < n; i++) {
        x.set(i, 0, 3.0 * (i + 1.0));
        for (size_t j = 0; j <= i; j++) {
            A.set(i, j, 4.0 * i + 7.0 * j + 1.0);
        }
    }
    double correct_val = 855904.5f;
    double val = A.quad(x);
    _ASSERT_EQ(correct_val, val);

    Matrix q(n, 1);
    for (size_t i = 0; i < n; i++) {
        q.set(i, 0, -5.0 * i - 1.0);
    }

    val = A.quad(x, q);
    correct_val = 850789.5;
    _ASSERT_EQ(correct_val, val);
}

void TestMatrix::testLeftTransposeMultiply() {
    size_t n = 10;
    size_t m = 5;
    size_t k = 8;

    Matrix A = MatrixFactory::MakeRandomMatrix(m, n, 00., 1.0, Matrix::MATRIX_DENSE);
    Matrix B = MatrixFactory::MakeRandomMatrix(m, k, 00., 1.0, Matrix::MATRIX_DENSE);

    A.transpose();
    _ASSERT_EQ(n, A.getNrows());
    _ASSERT_EQ(m, A.getNcols());

    Matrix C = A * B;

    _ASSERT_EQ(n, C.getNrows());
    _ASSERT_EQ(k, C.getNcols());

    A.transpose();
    Matrix AT(n, m, Matrix::MATRIX_DENSE);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            AT.set(i, j, A.get(j, i));
        }
    }
    Matrix C0 = AT * B;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            _ASSERT_EQ(C0.get(i, j), C.get(i, j));
        }
    }

}

void TestMatrix::testRightTransposeMultiply() {
    const size_t n = 10;
    const size_t m = 5;
    const size_t k = 8;

    Matrix A = MatrixFactory::MakeRandomMatrix(n, m, 00., 1.0, Matrix::MATRIX_DENSE);
    Matrix B = MatrixFactory::MakeRandomMatrix(k, m, 00., 1.0, Matrix::MATRIX_DENSE);

    B.transpose();
    Matrix C = A*B;

    B.transpose();
    Matrix BT(m, k);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < k; j++) {
            BT.set(i, j, B.get(j, i));
        }
    }

    Matrix C_correct = A*BT;
    _ASSERT_EQ(C_correct.getNrows(), C.getNrows());
    _ASSERT_EQ(C_correct.getNcols(), C.getNcols());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < k; j++) {
            _ASSERT_EQ(C_correct.get(i, j), C.get(i, j));
        }
    }
}

void TestMatrix::testLeftSymmetricMultiply() {
    const size_t n = 5;
    Matrix S = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
    Matrix A = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix C = S*A;

    Matrix S2D(n, n, Matrix::MATRIX_DENSE);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            S2D.set(i, j, S.get(i, j));
        }
    }

    Matrix C_correct = S2D*A;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(C_correct.get(i, j), C.get(i, j));
        }
    }

    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    C = S*x;

    C_correct = S2D * x;

    for (size_t i = 0; i < n; i++) {
        _ASSERT_NUM_EQ(C_correct.get(i, 0), C.get(i, 0), 1e-4);
    }
}

void TestMatrix::testSparseGetSet() {
    size_t n = 5;
    size_t m = 10;
    const size_t max_nnz = 3;
    Matrix M = MatrixFactory::MakeSparse(n, m, max_nnz, Matrix::SPARSE_UNSYMMETRIC);

    double r[3] = {4.576, 3.645, 1.092};
    M.set(0, 0, r[0]);
    M.set(0, 1, r[1]);
    M.set(1, 1, r[2]);

    _ASSERT_EQ(r[0], M.get(0, 0));
    _ASSERT_EQ(r[1], M.get(0, 1));
    _ASSERT_EQ(r[2], M.get(1, 1));
    _ASSERT_EQ(n, M.getNrows());
    _ASSERT_EQ(m, M.getNcols());
    _ASSERT_EQ(Matrix::MATRIX_SPARSE, M.getType());
    _ASSERT_EQ(0, M.cholmod_handle()->status);

    _ASSERT_OK(Matrix::destroy_handle());

    Matrix S = MatrixFactory::MakeSparse(n, m, max_nnz, Matrix::SPARSE_UNSYMMETRIC);
    for (size_t i = 0; i <= 20; i++) {
        _ASSERT_OK(S.set(0, 0, 0.5 * i));
    }
    _ASSERT_EQ(10.0, S.get(0, 0));

    _ASSERT_OK(S.set(1, 0, 1.0));
    _ASSERT_EQ(1.0, S.get(1, 0));

    _ASSERT_OK(S.set(2, 2, 5.0));
    _ASSERT_EQ(5.0, S.get(2, 2));



}

void TestMatrix::test_MSD() {

    size_t n = 3;
    size_t m = 3;
    size_t max_nnz = 4;

    Matrix A = MatrixFactory::MakeSparse(n, m, max_nnz, Matrix::SPARSE_UNSYMMETRIC);
    A.set(0, 0, 4.0);
    A.set(1, 0, 1.0); // A is declared as SPARSE_SYMMETRIC - no need to define A(0,1).
    A.set(1, 1, 5.0);
    A.set(2, 2, 10.0);


    Matrix B(m, n);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            B.set(i, j, 3.0 * i + 4.23 * j + 1.10);
        }
    }


    _ASSERT_EQ(m, B.getNrows());
    _ASSERT_EQ(n, B.getNcols());

    Matrix C;
    _ASSERT_OK(C = A * B);

    double correctData[9] = {
        4.4000, 21.6000, 71.0000,
        21.3200, 46.9800, 113.3000,
        38.2400, 72.3600, 155.6000
    };

    Matrix C_correct(n, n, correctData, Matrix::MATRIX_DENSE);


    _ASSERT_EQ(C_correct, C);
    _ASSERT_EQ(0, Matrix::cholmod_handle()->status);

}

void TestMatrix::test_MSDT() {
    const size_t n = 13;
    const size_t m = 25;
    const size_t nnz = 18;
    const size_t rep = 50;

    for (size_t k = 0; k < rep; k++) {
        Matrix S = MatrixFactory::MakeRandomSparse(m, n, nnz, 0.2, 2.0);
        Matrix D = MatrixFactory::MakeRandomMatrix(m, n, 0.5, 3.0, Matrix::MATRIX_DENSE);

        Matrix S2D(m, n);
        for (size_t i = 0; i < S.getNrows(); i++) {
            for (size_t j = 0; j < S.getNcols(); j++) {
                S2D.set(i, j, S.get(i, j));
            }
        }

        D.transpose();

        Matrix R = S*D; /* R = S * D' */
        Matrix correct = S2D*D;

        const double tol = 1e-8;
        for (size_t i = 0; i < R.getNrows(); i++) {
            for (size_t j = 0; j < R.getNcols(); j++) {
                _ASSERT_NUM_EQ(correct.get(i, j), R.get(i, j), tol);
            }
        }

    }
}

void TestMatrix::test_MSTDT() {
    const size_t n = 5;
    const size_t m = 8;
    const size_t nnz = 23;
    Matrix S = MatrixFactory::MakeRandomSparse(m, n, nnz, 0.2, 2.0);
    Matrix D = MatrixFactory::MakeRandomMatrix(n, m, 0.5, 3.0, Matrix::MATRIX_DENSE);

    Matrix DS = D*S;
    DS.transpose();

    S.transpose(); // n-m
    D.transpose(); // m-n

    Matrix R = S * D;

    /*
     * S'*D' = (D*S)'
     */
    _ASSERT_EQ(DS.getNrows(), R.getNrows());
    _ASSERT_EQ(DS.getNcols(), R.getNcols());
    for (size_t i = 0; i < R.getNrows(); i++) {
        for (size_t j = 0; j < R.getNcols(); j++) {
            _ASSERT_NUM_EQ(DS.get(i, j), R.get(i, j), 1e-8);
        }
    }

}

void TestMatrix::test_MDS() {
    const size_t n = 5;
    const size_t m = 7;

    Matrix D = MatrixFactory::MakeRandomMatrix(n, m, 0.5, 3.0, Matrix::MATRIX_DENSE);
    Matrix S = MatrixFactory::MakeRandomSparse(m, m, static_cast<size_t> (2.5 * (n + m)), 0.2, 2.0);

    Matrix R = D*S;

    D.transpose();
    S.transpose();

    Matrix C = S * D;

    C.transpose();

    _ASSERT_EQ(C.getNrows(), R.getNrows());
    _ASSERT_EQ(C.getNrows(), n);
    _ASSERT_EQ(C.getNcols(), R.getNcols());
    _ASSERT_EQ(C.getNcols(), m);
    for (size_t i = 0; i < C.getNrows(); i++) {
        for (size_t j = 0; j < C.getNcols(); j++) {
            _ASSERT_NUM_EQ(C.get(i, j), R.get(i, j), 1e-8);
        }
    }

}

void TestMatrix::test_MSS() {
    const size_t n = 10;
    const size_t m = 12;
    const size_t k = 9;
    const size_t nnz_L = 50;
    const size_t nnz_R = 40;

    Matrix L = MatrixFactory::MakeRandomSparse(n, m, nnz_L, 1.0, 2.0);
    Matrix R = MatrixFactory::MakeRandomSparse(m, k, nnz_R, 1.0, 2.0);

    Matrix Y;
    _ASSERT_OK(Y = L * R);
    _ASSERT_EQ(n, Y.getNrows());
    _ASSERT_EQ(k, Y.getNcols());
    _ASSERT_EQ(Matrix::MATRIX_SPARSE, Y.getType());

    Matrix Ld(n, m);
    Matrix Rd(m, k);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            Ld.set(i, j, L.get(i, j));
        }
    }

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < k; j++) {
            Rd.set(i, j, R.get(i, j));
        }
    }

    Matrix Z = Ld*Rd;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < k; j++) {
            _ASSERT_NUM_EQ(Z.get(i, j), Y.get(i, j), 1e-6);
        }
    }
}

void TestMatrix::testSparseAddDense() {
    const size_t n = 5;
    const size_t m = 7;
    const size_t nnz = 6;
    Matrix A = MatrixFactory::MakeSparse(n, m, nnz, Matrix::SPARSE_UNSYMMETRIC);
    A.set(0, 0, 0.5);
    A.set(0, 1, 0.2);
    A.set(1, 0, 1.2);
    A.set(1, 1, 3.6);
    A.set(4, 5, 6.2);
    A.set(3, 6, 9.9);

    Matrix B(n, m, Matrix::MATRIX_DENSE);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            B.set(i, j, -i - 3 * j - 1.5f);
        }
    }

    Matrix A_init(A);
    A += B;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            _ASSERT_EQ(A.get(i, j), A_init.get(i, j) + B.get(i, j));
        }
    }

}

void TestMatrix::testSparseAddSparse() {
    FILE *fp_A, *fp_B;
    fp_A = fopen("matrices/sparse2.mx", "r");
    fp_B = fopen("matrices/sparse3.mx", "r");

    CPPUNIT_ASSERT_MESSAGE("Can't open sparse2.mx", fp_A != NULL);
    CPPUNIT_ASSERT_MESSAGE("Can't open sparse3.mx", fp_B != NULL);

    Matrix A = MatrixFactory::ReadSparse(fp_A);
    Matrix B = MatrixFactory::ReadSparse(fp_B);

    /* close file handlers */
    _ASSERT_EQ(0, fclose(fp_A));
    _ASSERT_EQ(0, fclose(fp_B));

    Matrix A0 = A;
    _ASSERT_OK(A += B);

    for (size_t i = 0; i < A.getNrows(); i++) {
        for (size_t j = 0; j < A.getNcols(); j++) {
            _ASSERT_EQ(A0.get(i, j) + B.get(i, j), A.get(i, j));
        }
    }

}

void TestMatrix::testSparseAddSparse2() {
    const size_t n = 100;
    const size_t m = 250;
    const size_t nnz = 200;
    Matrix A = MatrixFactory::MakeRandomSparse(n, m, nnz, -1.0, 2.0);
    Matrix B = MatrixFactory::MakeRandomSparse(n, m, nnz, -1.0, 2.0);
    Matrix C;
    _ASSERT_OK(C = A + B);
    _ASSERT_EQ(Matrix::MATRIX_SPARSE, C.getType());
    _ASSERT_EQ(n, C.getNrows());
    _ASSERT_EQ(m, C.getNcols());

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            _ASSERT_EQ(A.get(i, j) + B.get(i, j), C.get(i, j));
        }
    }

    _ASSERT_OK(C += A);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            _ASSERT_EQ(2 * A.get(i, j) + B.get(i, j), C.get(i, j));
        }
    }

    _ASSERT_OK(Matrix::destroy_handle());
}

void TestMatrix::testSparseQuad() {
    size_t n = 60;
    size_t nnz = std::floor(1.2f * n);
    size_t tests = 20;

    Matrix *A = new Matrix();
    Matrix *x = new Matrix();

    for (size_t k = 0; k < tests; k++) {
        _ASSERT_OK(*A = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.0, 10.0));
        _ASSERT_OK(*x = MatrixFactory::MakeRandomMatrix(n, 1, 1.0, 2.0, Matrix::MATRIX_DENSE));

        double r = A->quad(*x);

        double r_exp = 0.0;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                r_exp += x->get(i, 0) * x->get(j, 0) * A->get(i, j);
            }
        }
        r *= 2;
        double tol = 1e-6;
        _ASSERT(std::abs(r) > tol);
        _ASSERT_NUM_EQ(r_exp, r, tol);
        _ASSERT_EQ(0, Matrix::cholmod_handle()->status);
    }

    delete A;
    delete x;

    _ASSERT_EQ(0, Matrix::cholmod_handle()->status);
    _ASSERT_OK(Matrix::destroy_handle());
}

void TestMatrix::testSparseQuadSparseX() {
    const double tol = 1e-6;
    const size_t runs = 10;
    for (size_t p = 0; p < runs; p++) {
        for (size_t n = 2; n < 70; n += 3) {
            size_t nnz_A = static_cast<size_t> (std::ceil(1.2 * n));
            size_t nnz_x = std::max(static_cast<size_t> (1), static_cast<size_t> (std::ceil(0.75 * n)));
            Matrix A = MatrixFactory::MakeRandomSparse(n, n, nnz_A, 0.0, 10.0);
            Matrix x = MatrixFactory::MakeRandomSparse(n, 1, nnz_x, 0.0, 10.0);
            double r, r_exp = 0.0;
            _ASSERT_OK(r = A.quad(x));
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    r_exp += x.get(i, 0) * x.get(j, 0) * A.get(i, j) / 2;
                }
            }
            _ASSERT_NUM_EQ(r_exp, r, tol);
        }
    }
}

void TestMatrix::testSparseQuad_q() {
    size_t n = 6;
    size_t nnz_A = 6;
    size_t nnz_x = 4;
    size_t nnz_q = 1;

    Matrix A = MatrixFactory::MakeSparse(n, n, nnz_A, Matrix::SPARSE_UNSYMMETRIC);
    A.set(0, 0, 4);
    A.set(0, 1, 5);
    A.set(1, 0, 8);
    A.set(1, 2, 10);
    A.set(5, 3, 100);

    Matrix x = MatrixFactory::MakeSparse(n, 1, nnz_x, Matrix::SPARSE_UNSYMMETRIC);
    x.set(0, 0, 1);
    x.set(1, 0, 2);
    x.set(5, 0, 9);
    x.set(4, 0, 3);

    Matrix q = MatrixFactory::MakeSparse(n, 1, nnz_q, Matrix::SPARSE_UNSYMMETRIC);
    q.set(4, 0, 1);

    double r = A.quad(x, q);
    const double r_exp = 18.0;
    const double tol = 1e-10;

    _ASSERT_NUM_EQ(r_exp, r, tol);
    _ASSERT_EQ(0, Matrix::cholmod_handle()->status);
    _ASSERT_OK(Matrix::destroy_handle());

}

void TestMatrix::testSparseDotProd() {
    const size_t n = 10;
    Matrix x = MatrixFactory::MakeSparse(n, 1, 0, Matrix::SPARSE_UNSYMMETRIC);
    Matrix result;
    _ASSERT_OK(result = x * x);
    _ASSERT_NOT(result.isEmpty());
    _ASSERT_EQ(static_cast<size_t> (1), result.getNcols());
    _ASSERT_EQ(static_cast<size_t> (1), result.getNrows());
    _ASSERT_EQ(0.0, result.get(0, 0));

    x = MatrixFactory::MakeSparse(n, 1, 1, Matrix::SPARSE_UNSYMMETRIC);
    x.set(0, 0, 2.0);
    result = x*x;
    _ASSERT_EQ(4.0, result.get(0, 0));
}

/* Tests R = DENSE + (?) */

void TestMatrix::test_ADD2() {
    size_t n = 80;
    size_t m = 5;
    const double tol = 1e-10;
    Matrix D1 = MatrixFactory::MakeRandomMatrix(n, m, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix D2 = MatrixFactory::MakeRandomMatrix(n, m, -5.0, 6.0, Matrix::MATRIX_DENSE);

    Matrix R;
    _ASSERT_OK(R = D1 + D2);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            _ASSERT_NUM_EQ(D1.get(i, j) + D2.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ADH() {
    size_t n = 100;
    const double tol = 1e-10;
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix H = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_SYMMETRIC);

    Matrix R;
    _ASSERT_OK(R = D + H);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(D.get(i, j) + H.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ADL() {
    size_t n = 80;
    const double tol = 1e-10;
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix L = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_LOWERTR);

    Matrix R;
    _ASSERT_OK(R = D + L);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(D.get(i, j) + L.get(i, j), R.get(i, j), tol);
            if (i < j) {
                _ASSERT_NUM_EQ(R.get(i, j), D.get(i, j), tol);
            }
        }
    }


}

void TestMatrix::test_ADS() {
    size_t n = 50;
    size_t m = 120;
    size_t nnz = 200;
    const double tol = 1e-10;
    Matrix D = MatrixFactory::MakeRandomMatrix(n, m, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix S = MatrixFactory::MakeRandomSparse(n, m, nnz, 0.0, 1.0);

    Matrix R;
    _ASSERT_OK(R = D + S);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            _ASSERT_NUM_EQ(D.get(i, j) + S.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ADW() {
    //std::cerr << "test_ADW:::EMPTY TEST!!!!\n";
}

void TestMatrix::test_ADX() {
    size_t n = 100;
    const double tol = 1e-10;
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix X = MatrixFactory::MakeRandomMatrix(n, n, 1.0, 10.0, Matrix::MATRIX_DIAGONAL);
    Matrix R;
    _ASSERT_OK(R = D + X);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(D.get(i, j) + X.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_EH() {
    size_t n = 50;
    Matrix H = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_SYMMETRIC);
    Matrix H2(H);
    Matrix H3 = H2;
    _ASSERT(H2.getNrows() == n);
    _ASSERT(H2.getNcols() == n);
    _ASSERT(H2.getType() == Matrix::MATRIX_SYMMETRIC);
    _ASSERT(H3.getNrows() == n);
    _ASSERT(H3.getNcols() == n);
    _ASSERT(H3.getType() == Matrix::MATRIX_SYMMETRIC);
    _ASSERT(H2.isSymmetric());
    _ASSERT(H3.isSymmetric());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(H.get(i, j), H2.get(i, j));
            _ASSERT_EQ(H.get(i, j), H3.get(i, j));
        }
    }
}

void TestMatrix::test_EL() {
    size_t n = 50;
    Matrix L = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_LOWERTR);
    Matrix L2(L);
    Matrix L3 = L2;
    _ASSERT(L2.getNrows() == n);
    _ASSERT(L2.getNcols() == n);
    _ASSERT(L2.getType() == Matrix::MATRIX_LOWERTR);
    _ASSERT(L3.getNrows() == n);
    _ASSERT(L3.getNcols() == n);
    _ASSERT(L3.getType() == Matrix::MATRIX_LOWERTR);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(L.get(i, j), L2.get(i, j));
            _ASSERT_EQ(L.get(i, j), L3.get(i, j));
        }
    }
}

void TestMatrix::test_EX() {
    size_t n = 50;
    Matrix X = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DIAGONAL);
    Matrix X2(X);
    Matrix X3 = X2;
    _ASSERT(X2.getNrows() == n);
    _ASSERT(X2.getNcols() == n);
    _ASSERT(X2.getType() == Matrix::MATRIX_DIAGONAL);
    _ASSERT(X3.getNrows() == n);
    _ASSERT(X3.getNcols() == n);
    _ASSERT(X3.getType() == Matrix::MATRIX_DIAGONAL);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(X.get(i, j), X2.get(i, j));
            _ASSERT_EQ(X.get(i, j), X3.get(i, j));
            if (i != j) {
                _ASSERT_EQ(0.0, X.get(i, j));
            }
        }
    }
}

void TestMatrix::test_EHT() {
    size_t n = 50;
    Matrix H = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_SYMMETRIC);
    Matrix H2(H);
    H2.transpose();
    Matrix H3 = H2;
    H3.transpose();
    _ASSERT(H2.getNrows() == n);
    _ASSERT(H2.getNcols() == n);
    _ASSERT(H2.getType() == Matrix::MATRIX_SYMMETRIC);
    _ASSERT(H3.getNrows() == n);
    _ASSERT(H3.getNcols() == n);
    _ASSERT(H3.getType() == Matrix::MATRIX_SYMMETRIC);
    _ASSERT(H2.isSymmetric());
    _ASSERT(H3.isSymmetric());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(H.get(i, j), H2.get(i, j));
            _ASSERT_EQ(H.get(i, j), H2.get(j, i));
            _ASSERT_EQ(H.get(i, j), H3.get(i, j));
            _ASSERT_EQ(H.get(i, j), H3.get(j, i));
        }
    }
}

void TestMatrix::test_EDT() {
    size_t n = 100;
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix D1 = D;
    Matrix D2(D);
    D1.transpose();
    D2.transpose();
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_EQ(D.get(i, j), D1.get(j, i));
            _ASSERT_EQ(D.get(i, j), D2.get(j, i));
        }
    }
}

void TestMatrix::test_ADDT() {
    size_t n = 10;
    const double tol = 1e-10;
    Matrix D1 = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix D2 = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix D2t(D2); // D2t = D2
    D2t.transpose(); // D2t = D2t'

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(D2.get(i, j), D2t.get(j, i), tol);
        }
    }

    Matrix R;
    _ASSERT_OK(R = D1 + D2t);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(D1.get(i, j) + D2.get(j, i), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ADHT() {
    size_t n = 10;
    const double tol = 1e-10;
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix H = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_SYMMETRIC);
    Matrix Ht(H); // D2t = D2
    Ht.transpose(); // D2t = D2t'

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(H.get(i, j), Ht.get(j, i), tol);
        }
    }

    Matrix R;
    _ASSERT_OK(R = D + Ht);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(D.get(i, j) + H.get(j, i), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ADLT() {
    size_t n = 10;
    const double tol = 1e-10;
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix L = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_LOWERTR);
    Matrix Lt(L); // D2t = D2
    Lt.transpose(); // D2t = D2t'

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(L.get(i, j), Lt.get(j, i), tol);
        }
    }

    Matrix R;
    _ASSERT_OK(R = D + Lt);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(D.get(i, j) + L.get(j, i), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ADST() {
    size_t n = 50;
    size_t nnz = 200;
    const double tol = 1e-10;
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix S = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.0, 10.0);
    Matrix St = S; // D2t = D2
    St.transpose(); // D2t = D2t'

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(S.get(i, j), St.get(j, i), tol);
        }
    }

    Matrix R;
    _ASSERT_OK(R = D + St);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(D.get(i, j) + S.get(j, i), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ADWT() {

}

void TestMatrix::test_ADXT() {
    size_t n = 10;
    const double tol = 1e-10;
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DENSE);
    Matrix X = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DIAGONAL);
    Matrix Xt(X); // D2t = D2
    Xt.transpose(); // D2t = D2t'

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(X.get(i, j), Xt.get(j, i), tol);
        }
    }

    Matrix R;
    _ASSERT_OK(R = D + Xt);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(D.get(i, j) + X.get(j, i), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_AHH() {
    size_t n = 100;
    const double tol = 1e-10;
    Matrix H1 = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_SYMMETRIC);
    Matrix H2 = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_SYMMETRIC);

    Matrix R;
    _ASSERT_OK(R = H1 + H2);
    _ASSERT_EQ(Matrix::MATRIX_SYMMETRIC, R.getType());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(H1.get(i, j) + H2.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_AHX() {
    size_t n = 100;
    const double tol = 1e-10;
    Matrix H = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_SYMMETRIC);
    Matrix X = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DIAGONAL);

    Matrix R;
    _ASSERT_EQ(Matrix::MATRIX_SYMMETRIC, H.getType());
    _ASSERT_EQ(Matrix::MATRIX_DIAGONAL, X.getType());
    _ASSERT_OK(R = H + X);
    _ASSERT_EQ(n, R.getNrows());
    _ASSERT_EQ(n, R.getNcols());
    _ASSERT(R.isSymmetric());
    _ASSERT(H.isSymmetric());
    _ASSERT(X.isSymmetric());
    _ASSERT_EQ(Matrix::MATRIX_SYMMETRIC, R.getType());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(H.get(i, j) + X.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_AHD() {
    size_t n = 100;
    const double tol = 1e-10;
    Matrix H = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_SYMMETRIC);
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_DENSE);

    Matrix R;
    _ASSERT_EQ(Matrix::MATRIX_SYMMETRIC, H.getType());
    R = H + D;
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType()); /* H + D = D */
    _ASSERT_EQ(n, R.getNrows());
    _ASSERT_EQ(n, R.getNcols());
    _ASSERT_NOT(R.isSymmetric());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(H.get(i, j) + D.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_AHL() {
    size_t n = 100;
    const double tol = 1e-10;
    Matrix H = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_SYMMETRIC);
    Matrix L = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_LOWERTR);

    Matrix R;
    _ASSERT_EQ(Matrix::MATRIX_SYMMETRIC, H.getType());
    _ASSERT_OK(R = H + L);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType()); /* H + L = D */
    _ASSERT_EQ(n, R.getNrows());
    _ASSERT_EQ(n, R.getNcols());
    _ASSERT_NOT(R.isSymmetric());
    _ASSERT(H.isSymmetric());
    _ASSERT_NOT(L.isSymmetric());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(H.get(i, j) + L.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_AHS() {
    size_t n = 40;
    size_t nnz = 90;
    const double tol = 1e-10;
    Matrix H = MatrixFactory::MakeRandomMatrix(n, n, -5.0, 6.0, Matrix::MATRIX_SYMMETRIC);
    Matrix S = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.0, 1.0);

    Matrix R;
    _ASSERT_EQ(Matrix::MATRIX_SYMMETRIC, H.getType());
    _ASSERT_OK(R = H + S);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType()); /* H + S = D */
    _ASSERT_EQ(n, R.getNrows());
    _ASSERT_EQ(n, R.getNcols());
    _ASSERT_NOT(R.isSymmetric());
    _ASSERT(H.isSymmetric());
    _ASSERT_NOT(S.isSymmetric());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(H.get(i, j) + S.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ASD() { /* Sparse + Diagonal */
    size_t n = 120;
    size_t nnz = 60;
    const double tol = 1e-10;
    Matrix S = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.0, 1.0);
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_DENSE);

    Matrix R;
    _ASSERT_OK(R = S + D);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());
    _ASSERT_EQ(n, R.getNrows());
    _ASSERT_EQ(n, R.getNcols());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(S.get(i, j) + D.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ASH() {
    size_t n = 120;
    size_t nnz = 60;
    const double tol = 1e-10;
    Matrix S = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.0, 1.0);
    Matrix H = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_SYMMETRIC);

    Matrix R;
    _ASSERT_OK(R = S + H);
    _ASSERT_EQ(Matrix::MATRIX_DENSE, R.getType());
    _ASSERT_EQ(n, R.getNrows());
    _ASSERT_EQ(n, R.getNcols());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(S.get(i, j) + H.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ASX() {
    size_t n = 120;
    size_t nnz = 60;
    const double tol = 1e-10;
    Matrix S = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.0, 1.0);
    Matrix X = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_DIAGONAL);

    Matrix R;
    _ASSERT_OK(R = S + X); // SPARSE + DIAGONAL = SPARSE
    _ASSERT_EQ(Matrix::MATRIX_SPARSE, R.getType());
    _ASSERT_EQ(n, R.getNrows());
    _ASSERT_EQ(n, R.getNcols());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(S.get(i, j) + X.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ASL() {
    size_t n = 20;
    size_t nnz = 60;
    const double tol = 1e-10;
    Matrix S = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.0, 1.0);
    Matrix L = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_LOWERTR);

    Matrix R;
    _ASSERT_OK(R = S + L); // SPARSE + LOWER TRIANGULAR = SPARSE
    _ASSERT_EQ(Matrix::MATRIX_SPARSE, R.getType());
    _ASSERT_EQ(n, R.getNrows());
    _ASSERT_EQ(n, R.getNcols());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(S.get(i, j) + L.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_AXX() {
    size_t n = 120;
    const double tol = 1e-10;
    Matrix X1 = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_DIAGONAL);
    Matrix X2 = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_DIAGONAL);

    Matrix R;
    _ASSERT_OK(R = X1 + X2); // DIAGONAL + DIAGONAL = DIAGONAL
    _ASSERT_EQ(Matrix::MATRIX_DIAGONAL, R.getType());
    _ASSERT_EQ(n, R.getNrows());
    _ASSERT_EQ(n, R.getNcols());
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(X1.get(i, j) + X2.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ALL() {
    size_t n = 120;
    const double tol = 1e-10;
    Matrix L1 = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_LOWERTR);
    Matrix L2 = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_LOWERTR);
    Matrix L;
    L = L1 + L2;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(L1.get(i, j) + L2.get(i, j), L.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ALX() {
    size_t n = 120;
    const double tol = 1e-10;
    Matrix L = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_LOWERTR);
    Matrix X = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_DIAGONAL);
    Matrix R;
    R = L + X;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(L.get(i, j) + X.get(i, j), R.get(i, j), tol);
        }
    }
}

void TestMatrix::test_CD() {
    size_t n = 120;
    const double tol = 1e-10;
    double alpha = 1.55;
    Matrix D = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix DD = D;
    D *= alpha;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(alpha * DD.get(i, j), D.get(i, j), tol);
        }
    }

}

void TestMatrix::test_CH() {
    size_t n = 120;
    const double tol = 1e-10;
    double alpha = -2.3456;
    Matrix H = MatrixFactory::MakeRandomMatrix(n, n, 2.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix HH = H;
    H *= alpha;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(alpha * HH.get(i, j), H.get(i, j), tol);
        }
    }
}

void TestMatrix::test_CS() {
    size_t n = 200;
    size_t nnz = 60;
    const double tol = 1e-10;
    double alpha = -2.3456;
    Matrix S = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.0, 1.0);
    Matrix SS = S;
    S *= alpha;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(alpha * SS.get(i, j), S.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ASST() {
    const size_t n = 10;
    const size_t m = n + 2;
    const size_t nnz = 2 * n + 10;
    const double tol = 1e-10;

    Matrix A = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.0, 1.0);
    Matrix At(A);
    At.transpose();

    A += At;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(A.get(i, j), A.get(j, i), tol);
        }
    }

    // and some more testing:
    const size_t repetitions = 20;
    Matrix B, C, Corig, X;
    for (size_t k = 0; k < repetitions; k++) {
        B = MatrixFactory::MakeRandomSparse(n, m, nnz, 0.0, 1.0);
        C = MatrixFactory::MakeRandomSparse(m, n, nnz, 0.0, 1.0);
        Corig = C;
        C.transpose(); // C is now n-by-m
        _ASSERT_EQ(n, C.getNrows());
        _ASSERT_EQ(m, C.getNcols());

        X = B + C;

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < m; j++) {
                _ASSERT_NUM_EQ(X.get(i, j), B.get(i, j) + Corig.get(j, i), tol);
            }
        }
    }

}

void TestMatrix::test_EST() {
    const size_t n = 100;
    const size_t m = 80;
    const size_t nnz = 300;
    Matrix S = MatrixFactory::MakeRandomSparse(n, m, nnz, 0.0, 1.0);
    Matrix S1 = S;
    Matrix S2(S);
    S1.transpose();
    S2.transpose();
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            _ASSERT_EQ(S.get(i, j), S1.get(j, i));
            _ASSERT_EQ(S.get(i, j), S2.get(j, i));
        }
    }
}

void TestMatrix::test_ES() {
    const size_t n = 100;
    const size_t m = 80;
    const size_t nnz = 300;
    Matrix S = MatrixFactory::MakeRandomSparse(n, m, nnz, 0.0, 1.0);
    Matrix S1 = S;
    Matrix S2(S);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            _ASSERT_EQ(S.get(i, j), S1.get(i, j));
            _ASSERT_EQ(S.get(i, j), S2.get(i, j));
        }
    }
}

void TestMatrix::test_ASS() {
    const double tol = 1e-8;
    const size_t n = 100;
    const size_t m = 80;
    const size_t nnz = 300;
    Matrix A = MatrixFactory::MakeRandomSparse(n, m, nnz, 0.0, 1.0);
    Matrix B = MatrixFactory::MakeRandomSparse(n, m, nnz, 0.0, 1.0);
    Matrix X = A + B;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            _ASSERT_NUM_EQ(A.get(i, j) + B.get(i, j), X.get(i, j), tol);
        }
    }
}

void TestMatrix::testSubmatrix() {
    const double tol = 1e-12;

    const size_t n = 14;
    const size_t m = 11;

    const size_t row_start = 1;
    const size_t row_end = 1;
    const size_t col_start = 8;
    const size_t col_end = 8;

    Matrix A = MatrixFactory::MakeRandomMatrix(m, n, 2.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix subA = A.submatrixCopy(row_start, row_end, col_start, col_end);
    _ASSERT_EQ(row_end - row_start + 1, subA.getNrows());
    _ASSERT_EQ(col_end - col_start + 1, subA.getNcols());

    for (size_t i = 0; i < subA.getNrows(); i++) {
        for (size_t j = 0; j < subA.getNcols(); j++) {
            _ASSERT_NUM_EQ(A.get(row_start + i, col_start + j), subA.get(i, j), tol);
        }
    }
}

void TestMatrix::testSubmatrixTranspose() {
    const size_t n = 15;
    const size_t m = 12;
    const double tol = 1e-12;
    Matrix A = MatrixFactory::MakeRandomMatrix(m, n, 2.0, 10.0, Matrix::MATRIX_DENSE);

    A.transpose();

    size_t row_start_tr;
    size_t col_start_tr;
    size_t rows_tr;
    size_t cols_tr;
    Matrix subA;

    for (rows_tr = 0; rows_tr <= 4; rows_tr++) {
        for (cols_tr = 0; cols_tr <= 4; cols_tr++) {
            for (row_start_tr = 0; row_start_tr < n - rows_tr - 1; row_start_tr++) {
                for (col_start_tr = 0; col_start_tr < m - cols_tr - 1; col_start_tr++) {
                    subA = A.submatrixCopy(row_start_tr, row_start_tr + rows_tr, col_start_tr, col_start_tr + cols_tr);
                    for (size_t i = 0; i < subA.getNrows(); i++) {
                        for (size_t j = 0; j < subA.getNcols(); j++) {
                            _ASSERT_NUM_EQ(A.get(row_start_tr + i, col_start_tr + j), subA.get(i, j), tol);
                        }
                    }
                }
            }
        }
    }

}

void _testSubmatrixMultiply(Matrix& A, Matrix& B) {
    Matrix Asub;
    Matrix Bsub;
    Matrix exact;
    Matrix result;


    for (size_t k = 0; k < 3; k++) {
        for (size_t ds = 0; ds < 3; ds++) {
            for (size_t s = 0; s < A.getNcols() - ds - 1; s++) {
                for (size_t dj = 0; dj < 6; dj++) {
                    for (size_t j = 0; j < B.getNcols() - dj - 1; j++) {
                        for (size_t di = 0; di < 2; di++) {
                            for (size_t i = 0; i < A.getNrows() - di - 1; i++) {
                                Asub = A.submatrixCopy(i, i + di, s, s + ds);
                                Bsub = B.submatrixCopy(k, k + ds, j, j + dj);
                                exact = Asub*Bsub;
                                result = A.multiplySubmatrix(B,
                                        i, i + di, s, s + ds,
                                        k, k + ds, j, j + dj);
                                _ASSERT_EQ(exact, result);
                                _ASSERT_EQ(di + 1, result.getNrows());
                                _ASSERT_EQ(dj + 1, result.getNcols());
                            }
                        }
                    }
                }
            }
        }
    }
}

void TestMatrix::testSubmatrixMultiply() {
    const size_t MA = 11;
    const size_t NA = 10;
    const size_t MB = 11;
    const size_t NB = 9;
    Matrix A = MatrixFactory::MakeRandomMatrix(MA, NA, 0.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix B = MatrixFactory::MakeRandomMatrix(MB, NB, 0.0, 2.0, Matrix::MATRIX_DENSE);
    _testSubmatrixMultiply(A, B);
}

void TestMatrix::testSubmatrixMultiplyTr() {
    const size_t MA = 7;
    const size_t NA = 9;
    const size_t MB = 8;
    const size_t NB = 10;
    Matrix A = MatrixFactory::MakeRandomMatrix(MA, NA, 0.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix B = MatrixFactory::MakeRandomMatrix(MB, NB, 0.0, 2.0, Matrix::MATRIX_DENSE);
    A.transpose();
    _testSubmatrixMultiply(A, B); // A is transposed
    A.transpose();
    B.transpose();
    _testSubmatrixMultiply(A, B); // B is transposed
    A.transpose();
    _testSubmatrixMultiply(A, B); // both A and B are transposed
}

void TestMatrix::testSubmatrixSparse() {
    const double tol = 1e-12;

    const size_t n = 14;
    const size_t m = 11;

    const size_t row_start = 1;
    const size_t row_end = 1;
    const size_t col_start = 8;
    const size_t col_end = 8;

    Matrix A = MatrixFactory::MakeRandomSparse(m, n, m + n, 0.0, 1.0);
    Matrix subA = A.submatrixCopy(row_start, row_end, col_start, col_end);
    _ASSERT_EQ(row_end - row_start + 1, subA.getNrows());
    _ASSERT_EQ(col_end - col_start + 1, subA.getNcols());

    for (size_t i = 0; i < subA.getNrows(); i++) {
        for (size_t j = 0; j < subA.getNcols(); j++) {
            _ASSERT_NUM_EQ(A.get(row_start + i, col_start + j), subA.get(i, j), tol);
        }
    }
}

void TestMatrix::test_ADTDT() {
    size_t n = 15;
    size_t m = 25;

    Matrix D1 = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix D2 = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix C1(D1);
    Matrix C2(D2);
    D1.transpose();
    D2.transpose();

    Matrix R = D1 + D2;
    Matrix C = C1 + C2;
    C.transpose();
    _ASSERT_EQ(R, C);

}

void TestMatrix::testGetSetTranspose() {
    size_t n = 5;
    size_t m = 6;
    Matrix A(n, m);
    A.set(1, 3, 10.4);

    A.transpose();
    _ASSERT_EQ(10.4, A.get(3, 1));

    A.set(4, 3, 0.5);
    _ASSERT_EQ(0.5, A.get(4, 3));
    _ASSERT_EQ(0.0, A.get(3, 4));
}

void TestMatrix::test_ASTDT() {
    size_t n = 5;
    size_t m = 6;

    Matrix D1 = MatrixFactory::MakeRandomSparse(n, m, 30, 0.0, 1.0);
    Matrix D2 = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0, Matrix::MATRIX_DENSE);

    D1.transpose();
    D2.transpose();

    Matrix F = D1 + D2;

    _ASSERT_EQ(m, F.getNrows());
    _ASSERT_EQ(n, F.getNcols());

    const double tol = 1e-8;
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            _ASSERT_NUM_EQ(D1.get(i, j) + D2.get(i, j), F.get(i, j), tol);
        }
    }
}



