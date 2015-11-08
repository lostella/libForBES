#include "ForBES.h"

int main(int argc, char** argv) {
    size_t n = 10;
    size_t m = 3;

    Matrix M = MatrixFactory::MakeRandomMatrix(m, n, 0.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix x = MatrixFactory::MakeRandomMatrix(n, 1, 0.0, 10.0, Matrix::MATRIX_DENSE);
    LinearOperator *T = new MatrixOperator(M);
    Matrix y;
    y = T->call(x);   // y = T(x)
    Matrix z = M*x;   // z = M*x

    // z and y are equal

    std::cout << "Result, y = \n" << y;
    std::cout << "\n\nz = \n" << z;
    delete T;
}

/* 
 * Typical output:
 * 
 * Result, y = 
 * Matrix 3x1
 * Type: Dense
 *  310.8
 *  374.4
 *  263.9
 *  
 * z = 
 * 
 * Matrix 3x1
 * Type: Dense
 *   310.8
 *   374.4
 *   263.9
 *
 */