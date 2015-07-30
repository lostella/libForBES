/* 
 * File:   main.cpp
 * Author: Chung
 *
 * Created on July 7, 2015, 7:47 PM
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>

#include "Function.h"
#include "Quadratic.h"

#include "cholmod.h"
#include "MatrixFactory.h"
#include "LDLFactorization.h"

#include <set>


using namespace std;

int main(int argc, char** argv) {

    const size_t n = 3;
    double a[n * n] = {
        10, 20, 30,
        20, 45, 80,
        30, 80, 171
    };

    Matrix A(n, n, a, Matrix::MATRIX_DENSE);

    int j = 0;
    for (size_t i = 0; i < n * n; i++) {
        std::cout << std::setw(8) << std::setprecision(4) << a[i] << " ";
        j++;
        if (j == 3) {
            std::cout << std::endl;
            j = 0;
        }
    }
    std::cout << std::endl;
    std::cout << std::endl;

    int ipiv[n];
    int status = LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'L', n, a, n, ipiv);


    for (size_t i = 0; i < n * n; i++) {
        std::cout << std::setw(8) << std::setprecision(4) << a[i] << " ";
        j++;
        if (j == 3) {
            std::cout << std::endl;
            j = 0;
        }
    }

    std::cout << std::endl;
    std::cout << std::endl;

    for (size_t i = 0; i < n; i++) {
        std::cout << ipiv[i] << " ";
    }

    double b[n] = {
        1.2, 3.4, -0.5
    };
    Matrix B(n, 1, b, Matrix::MATRIX_DENSE);

        int status2 = LAPACKE_dsytrs(LAPACK_COL_MAJOR, 'L', n, 1, a, n, ipiv, b, n);


    std::cout << std::endl;
    std::cout << std::endl;
    for (size_t i = 0; i < n; i++) {
        std::cout << "x[" << i << "] = " << b[i] << "\n";
    }
    Matrix X(n, 1, b, Matrix::MATRIX_DENSE);

    Matrix ERR = B - A*X;
    std::cout << ERR;
    
    FactoredSolver *ldlSolver = new LDLFactorization(A);
    ldlSolver->factorize();
    Matrix SOL;
    ldlSolver->solve(B, SOL);
    
    Matrix ERR2 = B - A*SOL;
    std::cout << ERR2;
    
    return (0);
}

