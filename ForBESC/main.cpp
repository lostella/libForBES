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
    const double tol = 1e-7;
    size_t n = 5;
    size_t nnz = 2 * n - 1;
    Matrix A = MatrixFactory::MakeSparseSymmetric(n, nnz);    
    Matrix b(n, 1);

    for (size_t i = 0; i < n; i++) {
        A.set(i, i, 10.0);
        b.set(i, 0, i + 1);
    }
    for (size_t i = 1; i < n; i++) { /* Set the LT part only */
        A.set(i, i - 1, 0.5);
    }


    Matrix x;
    CholeskyFactorization * solver = new CholeskyFactorization(A);
    solver -> factorize();
    solver -> solve(b, x);

    std::cout << b;
    
    Matrix ax = A * x;
    for (size_t i = 0; i < n; i++) {
        std::cout << ax.get(i,0) << ", " << (std::abs(ax.get(i, 0) - i - 1.0)) << "\n";
    }




    return (0);
}

