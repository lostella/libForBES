#include "ForBES.h"
#include <iostream>
#include <stdlib>

int main(int argc, char** argv) {
    
    const size_t n = 15;
    const size_t m = 5;

    // create a vector x: (15 x 1)
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 10.0, 1.0);
    // create a vector y: (15 x 1)
    Matrix y = MatrixFactory::MakeRandomMatrix(n, 1, -1.0, 1.0);

    // random matrix A: (5 x 5)
    Matrix A = MatrixFactory::MakeRandomMatrix(m, m, 0.0, 1.0);

    // reference x[3:7]
    Matrix x_sub = MatrixFactory::ShallowVector(x, m, 3);
    // reference y[8:12]
    Matrix y_sub = MatrixFactory::ShallowVector(y, m, 8);

    // perform y[8:12] = A * x[3:7]
    Matrix::add(y_sub, 2.0, x_sub, 1.0);
    
    std::cout << y;

    return EXIT_SUCCESS;
}
