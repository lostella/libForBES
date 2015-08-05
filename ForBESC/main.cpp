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
#include "CholeskyFactorization.h"

#include <set>


using namespace std;

int main(int argc, char** argv) {
    size_t n = 8;
    Matrix A = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
    Matrix At = A;
    At.transpose();
    A += At;
    A *= 0.5;
    for (int i = 0; i < n; i++) {
        A.set(i, i, A.get(i, i) + 50.0);
    }
    std::cout << A;

    Matrix b = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);

    CholeskyFactorization * cholFactorization = new CholeskyFactorization(A);
    int status = cholFactorization ->factorize();
    std::cout << "status = " << status << "\n";

    Matrix x;
    int status2 = cholFactorization ->solve(b, x);

    Matrix err = A * x - b;

    std::cout << err;


    delete cholFactorization;

    return (0);
}

