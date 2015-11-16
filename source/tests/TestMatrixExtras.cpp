/*
 * File:   TestMatrixExtras.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Nov 8, 2015, 4:33:46 PM
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

#include "TestMatrixExtras.h"
#include "MatrixFactory.h"
#include "ForBES.h"
#include <cmath>

CPPUNIT_TEST_SUITE_REGISTRATION(TestMatrixExtras);

TestMatrixExtras::TestMatrixExtras() {
}

TestMatrixExtras::~TestMatrixExtras() {
}

void TestMatrixExtras::setUp() {
}

void TestMatrixExtras::tearDown() {
    Matrix::destroy_handle();
}

void TestMatrixExtras::test_add_wrong_args() {
    Matrix A(5, 6);
    Matrix B(5, 7);
    _ASSERT_EXCEPTION(Matrix::add(A, 1.3, B, 0.8), std::invalid_argument);
}

void TestMatrixExtras::test_mult_wrong_args() {
    {
        Matrix A(5, 6);
        Matrix B(6, 7);
        Matrix C(5, 4);
        _ASSERT_EXCEPTION(Matrix::mult(C, 1.3, A, B, 0.8), std::invalid_argument);
    }
    {
        Matrix A(5, 6);
        Matrix B(3, 7);
        Matrix C(5, 7);
        _ASSERT_EXCEPTION(Matrix::mult(C, 1.3, A, B, 0.8), std::invalid_argument);
    }
}

/*****    ADDITION    ******/

void TestMatrixExtras::test_add_DD() {
    size_t n = 10;
    size_t m = 15;
    size_t repetitions = 20;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, m, 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(n, m, 2.0, 1.0);

        double alpha = 2.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
        double gamma = -0.5 + static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);

        Matrix A_copy(A);
        Matrix B_copy(B);

        // A = gamma*A + alpha*B
        int status = Matrix::add(A, alpha, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
        _ASSERT_EQ(B_copy, B);

        A_copy *= gamma;
        B_copy *= alpha;
        Matrix C = A_copy + B_copy;

        _ASSERT_EQ(C, A);
    }

}

void TestMatrixExtras::test_add_DS() {
    size_t n = 10;
    size_t m = 15;
    size_t nnz = 100;
    Matrix A = MatrixFactory::MakeRandomMatrix(n, m, 2.0, 1.0);
    Matrix B = MatrixFactory::MakeRandomSparse(n, m, nnz, 2.0, 1.0);

    double alpha = -1.0 + 2.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
    double gamma = -0.5 + static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);


    Matrix R(A);
    R *= gamma;
    Matrix T(B);
    T *= alpha;
    Matrix S = R + T;

    int status = Matrix::add(A, alpha, B, gamma); // A = gamma * A + alpha * B
    _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
    _ASSERT_EQ(S, A);

}

void TestMatrixExtras::test_add_DST() {
    for (size_t n = 11; n < 40; n += 5) {
        for (size_t m = 13; m < 30; m += 3) {
            size_t nnz = 100;
            Matrix A = MatrixFactory::MakeRandomMatrix(n, m, 2.0, 1.0);
            Matrix B = MatrixFactory::MakeRandomSparse(m, n, nnz, 2.0, 1.0);

            Matrix A_copy(A);
            Matrix B_copy(B);

            B.transpose();

            double alpha = -1.0 + 2.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
            double gamma = -0.5 + static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
            const double tol = 1e-8;

            int status = Matrix::add(A, alpha, B, gamma); // A = gamma * A + alpha * B
            _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

            _ASSERT_EQ(n, A.getNrows());
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < m; j++) {
                    _ASSERT_NUM_EQ(gamma * A_copy.get(i, j) + alpha * B_copy.get(j, i), A.get(i, j), tol);
                }
            }
        }
    }
}

void TestMatrixExtras::test_add_SS() {
    size_t n = 10;
    size_t m = 15;
    size_t nnz = 100;
    size_t repetitions = 300;

    Matrix A = MatrixFactory::MakeRandomSparse(n, m, nnz, 2.0, 1.0);
    for (size_t r = 0; r < repetitions; r++) {
        Matrix B = MatrixFactory::MakeRandomSparse(n, m, nnz, 2.0, 1.0);

        double alpha = 2.0 * static_cast<double> (std::rand()) / RAND_MAX;
        double gamma = -0.5 + static_cast<double> (std::rand()) / RAND_MAX;

        Matrix A_copy(A);
        Matrix B_copy(B);

        // A = gamma*A + alpha*B
        int status = Matrix::add(A, alpha, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
        _ASSERT_EQ(B_copy, B);

        A_copy *= gamma;
        B_copy *= alpha;
        Matrix C = A_copy + B_copy;

        _ASSERT_EQ(C, A);
    }
}

void TestMatrixExtras::test_add_DDT() {
    size_t n = 10;
    size_t m = 15;
    size_t repetitions = 20;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(m, n, 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(n, m, 2.0, 1.0);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        A.transpose();

        Matrix A_copy(A);
        Matrix B_copy(B);

        // A = gamma*A + alpha*B
        int status = Matrix::add(A, alpha, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
        _ASSERT_EQ(B_copy, B);

        _ASSERT_EQ(n, A.getNrows());
        _ASSERT_EQ(m, A.getNcols());

        A_copy *= gamma;
        B_copy *= alpha;
        Matrix C = A_copy + B_copy;

        _ASSERT_EQ(C, A);
    }
}

void TestMatrixExtras::test_add_DTD() {
    size_t n = 10;
    size_t m = 15;
    size_t repetitions = 20;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, m, 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(m, n, 2.0, 1.0);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        B.transpose();

        Matrix A_copy(A);
        Matrix B_copy(B);

        int status = Matrix::add(A, alpha, B, gamma); // A = gamma*A + alpha*B'
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
        _ASSERT_EQ(B_copy, B);

        _ASSERT_EQ(n, A.getNrows());
        _ASSERT_EQ(m, A.getNcols());

        A_copy *= gamma;
        B_copy *= alpha;
        Matrix C = A_copy + B_copy;

        _ASSERT_EQ(C, A);
    }
}

void TestMatrixExtras::test_add_DTDT() {
    size_t n = 10;
    size_t m = 15;
    size_t repetitions = 20;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, m, 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(n, m, 2.0, 1.0);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        Matrix A_copy(A);
        Matrix B_copy(B);

        A.transpose();
        B.transpose();

        // A := gamma*A + alpha*B
        int status = Matrix::add(A, alpha, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
        _ASSERT_EQ(m, A.getNrows());
        _ASSERT_EQ(n, A.getNcols());

        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                _ASSERT_EQ(gamma * A_copy.get(j, i) + alpha * B_copy.get(j, i), A.get(i, j));
            }
        }


    }
}

void TestMatrixExtras::test_add_SST() {
    size_t n = 10;
    size_t m = 15;
    size_t nnz = std::ceil(n * m / 2);
    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomSparse(n, m, nnz, 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomSparse(m, n, nnz, 2.0, 1.0);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        Matrix B_copy(B);
        B.transpose();
        Matrix A_copy(A);

        int status = Matrix::add(A, alpha, B, gamma); // A = gamma*A + alpha*B'
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        _ASSERT_EQ(n, A.getNrows());
        _ASSERT_EQ(m, A.getNcols());

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < m; j++) {
                _ASSERT_NUM_EQ(gamma * A_copy.get(i, j) + alpha * B_copy.get(j, i), A.get(i, j), 1e-8);
            }
        }

        Matrix S(A_copy);
        S *= gamma;
        Matrix T(B);
        T *= alpha;

        S += T;

        _ASSERT_EQ(S, A);

    }
}

void TestMatrixExtras::test_add_STS() {
    size_t n = 10;
    size_t m = 15;
    size_t nnz = 100;
    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix B = MatrixFactory::MakeRandomSparse(n, m, nnz, 2.0, 1.0);
        Matrix A = MatrixFactory::MakeRandomSparse(m, n, nnz, 2.0, 1.0);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);


        Matrix A_copy(A);
        A.transpose();

        int status = Matrix::add(A, alpha, B, gamma); // A' = gamma*A' + alpha*B
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        _ASSERT_EQ(n, A.getNrows());
        _ASSERT_EQ(m, A.getNcols());

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < m; j++) {
                _ASSERT_EQ(gamma * A_copy.get(j, i) + alpha * B.get(i, j), A.get(i, j));
            }
        }

    }
}

void TestMatrixExtras::test_add_STST() {
    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        size_t n = 10;
        size_t m = 15;
        size_t nnz = 70;
        Matrix A = MatrixFactory::MakeRandomSparse(n, m, nnz, 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomSparse(n, m, nnz, 2.0, 1.0);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);


        Matrix A_copy(A);
        Matrix B_copy(B);
        A.transpose();
        B.transpose();

        int status = Matrix::add(A, alpha, B, gamma); // A' = gamma*A' + alpha*B
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        _ASSERT_EQ(m, A.getNrows());
        _ASSERT_EQ(n, A.getNcols());
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                _ASSERT_EQ(gamma * A_copy.get(j, i) + alpha * B_copy.get(j, i), A.get(i, j));
            }
        }

    }
}

/*****    MULTIPLICATION ******/

void TestMatrixExtras::test_mult_DD() {
    size_t n = 8;
    size_t k = 6;
    size_t m = 5;
    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(k, m, 0.0, 1.0);
        Matrix C = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0);
        Matrix C_copy(C);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        int status = Matrix::mult(C, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        Matrix AB = A*B;
        AB *= alpha;
        C_copy *= gamma;
        C_copy += AB;

        _ASSERT_EQ(C_copy, C);
    }
}

void TestMatrixExtras::test_mult_DH() {
    size_t n = 8;
    size_t k = 6;
    size_t repetitions = 100;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(k, k, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
        Matrix C = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix C_copy(C);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        int status = Matrix::mult(C, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        Matrix AB = A*B;
        AB *= alpha;
        C_copy *= gamma;
        _ASSERT_NOT(C_copy == C);
        C_copy += AB;

        _ASSERT_EQ(C_copy, C);
    }
}

void TestMatrixExtras::test_mult_DDT() {
    size_t n = 8;
    size_t k = 6;
    size_t m = 5;
    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(m, k, 0.0, 1.0);
        Matrix C = MatrixFactory::MakeRandomMatrix(n, m, 0.0, 1.0);
        Matrix C_copy(C);

        B.transpose();

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        int status = Matrix::mult(C, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        Matrix AB = A*B;
        AB *= alpha;
        C_copy *= gamma;
        _ASSERT_NOT(C_copy == C);
        C_copy += AB;

        _ASSERT_EQ(C_copy, C);
    }
}

void TestMatrixExtras::test_mult_DX() {
    size_t n = 8;
    size_t k = 6;
    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(k, k, 0.0, 1.0, Matrix::MATRIX_DIAGONAL);
        Matrix C = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);
        Matrix C_copy(C);

        double alpha = 2.0 * static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);
        double gamma = -0.5 + static_cast<float> (std::rand()) / static_cast<float> (RAND_MAX);

        int status = Matrix::mult(C, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        Matrix AB = A*B;
        AB *= alpha;
        C_copy *= gamma;
        C_copy += AB;

        _ASSERT_EQ(C_copy, C);
    }
}

void TestMatrixExtras::test_mult_SS() {
    size_t n = 10;
    size_t k = 8;
    size_t m = 9;

    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {

        Matrix A = MatrixFactory::MakeRandomSparse(n, k, std::ceil(n * k / 2), 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomSparse(k, m, std::ceil(k * m / 2), 2.0, 1.0);
        Matrix C = MatrixFactory::MakeRandomSparse(n, m, std::ceil(n * m / 2), 2.0, 1.0);

        double alpha = -2.0 + 2.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
        double gamma = -2.0 + 3.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);


        Matrix aAB = A*B;
        aAB *= alpha; /*     aAB = a * A * B          */

        Matrix gC(C);
        gC *= gamma; /*     gC  = g * C              */

        Matrix R = aAB + gC;

        Matrix C_copy(C);
        int status = Matrix::mult(C_copy, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        _ASSERT_EQ(R, C_copy);
    }

}

void TestMatrixExtras::test_mult_SS2() {
    // testing with gamma = 0.0;
    size_t repetitions = 300;
    for (size_t r = 0; r < repetitions; r++) {
        size_t n = 10;
        size_t k = 8;
        size_t m = 9;
        Matrix A = MatrixFactory::MakeRandomSparse(n, k, std::ceil(n * k / 2), 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomSparse(k, m, std::ceil(k * m / 2), 2.0, 1.0);
        Matrix C = MatrixFactory::MakeRandomSparse(n, m, std::ceil(n * m / 2), 2.0, 1.0);

        double alpha = -2.0 + 2.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
        double gamma = 0.0;


        Matrix aAB = A*B;
        aAB *= alpha; /*     aAB = a * A * B          */


        Matrix C_copy(C);
        int status = Matrix::mult(C_copy, alpha, A, B, gamma); /*    C_copy = a * A * B        */
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);
        _ASSERT_EQ(aAB, C_copy);
    }
}

void TestMatrixExtras::test_mult_SS3() {
    // testing with C : DENSE

    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        size_t n = 10;
        size_t k = 8;
        size_t m = 9;
        Matrix A = MatrixFactory::MakeRandomSparse(n, k, std::ceil(n * k / 2), 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomSparse(k, m, std::ceil(m * k / 2), 2.0, 1.0);
        Matrix C = MatrixFactory::MakeRandomMatrix(n, m, 2.0, 1.0);

        double alpha = -2.0 + 2.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
        double gamma = -2.0 + 3.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);


        Matrix aAB = A*B;
        aAB *= alpha; /*     aAB = a * A * B          */

        Matrix gC(C);
        gC *= gamma; /*     gC  = g * C              */

        Matrix R = aAB + gC;

        Matrix C_copy(C);
        int status = Matrix::mult(C_copy, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        _ASSERT_EQ(R, C_copy);
    }
}

void TestMatrixExtras::test_mult_SST() {
    size_t n = 10;
    size_t k = 8;
    size_t m = 9;

    size_t repetitions = 300;

    for (size_t r = 0; r < repetitions; r++) {
        Matrix A = MatrixFactory::MakeRandomSparse(n, k, std::ceil(n * k / 2), 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomSparse(m, k, std::ceil(m * k / 2), 2.0, 1.0);
        Matrix C = MatrixFactory::MakeRandomSparse(n, m, std::ceil(n * m / 2), 2.0, 1.0);

        B.transpose();

        double alpha = -2.0 + 2.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
        double gamma = -2.0 + 3.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);

        Matrix C_orig(C);

        int status = Matrix::mult(C, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        C_orig *= gamma;
        Matrix AB = A * B;
        AB *= alpha;
        C_orig += AB;

        _ASSERT_EQ(C_orig, C);
    }
}

void TestMatrixExtras::test_mult_SX() {

    size_t repetitions = 300;
    for (size_t r = 0; r < repetitions; r++) {
        size_t n = 10;
        size_t k = 8;
        size_t nnz = std::ceil(n * k / 2);
        Matrix A = MatrixFactory::MakeRandomSparse(n, k, nnz, 2.0, 1.0);
        Matrix B = MatrixFactory::MakeRandomMatrix(k, k, 0.0, 1.0, Matrix::MATRIX_DIAGONAL);
        Matrix C = MatrixFactory::MakeRandomSparse(n, k, nnz, 2.0, 1.0);
        Matrix D = MatrixFactory::MakeRandomMatrix(n, k, 0.0, 1.0);

        double alpha = -2.0 + 2.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);
        double gamma = -2.0 + 3.0 * static_cast<double> (std::rand()) / static_cast<double> (RAND_MAX);

        Matrix D_orig(D);
        Matrix C_orig(C);

        int status = Matrix::mult(C, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        status = Matrix::mult(D, alpha, A, B, gamma);
        _ASSERT_EQ(ForBESUtils::STATUS_OK, status);

        C_orig *= gamma;
        D_orig *= gamma;
        Matrix AB = A * B;
        AB *= alpha;
        C_orig += AB;
        D_orig += AB;

        _ASSERT_EQ(C_orig, C);
        _ASSERT_EQ(D_orig, D);

    }
}

void TestMatrixExtras::test_mult_Hv() {
    size_t n = 10;
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0);
    Matrix A = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
    Matrix y = A*x;

    Matrix y2(n, 1);
    Matrix::mult(y2, 1.0, A, x, 0.0);

    _ASSERT_EQ(y2, y);
}
