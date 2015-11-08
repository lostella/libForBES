#include "ForBES.h"
#include <iostream>

int main(int argc, char** argv) {

    // define sizes
    size_t n = 5;
    size_t m = 7;

    // Create a dense matrix 5-by-7 with all entries set to 0.0:
    Matrix A(n, m, Matrix::MATRIX_DENSE);

    size_t nrows = A.getNrows();
    size_t ncols = A.getNcols();

    std::cout << "Matrix A is " << nrows << " x " << ncols << std::endl;

    // Set A(1,0) to a value:
    A.set(1, 0, 10.56);

    // retrieve the value of A(1,0):
    double value = A.get(1, 0);

    std::cout << "value = " << value << std::endl;

    // Make a random symmetric matrix:
    Matrix S = MatrixFactory::MakeRandomMatrix(n, n, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);

    // Make a random lower triangular matrix and random entries in [1, 3]:
    Matrix Ltri = MatrixFactory::MakeRandomMatrix(n, n, 1.0, 2.0, Matrix::MATRIX_LOWERTR);

    // Make a random sparse matrix with given number of non-zero entries:
    size_t nnz = 10;
    Matrix Y = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.01, 1.5);

    std::cout << "Y is a sparse matrix:\n";
    std::cout << Y;
    

    /*
     * Operations with matrices (addition, subtraction)
     */
    n = 10;
    nnz = 20;
    Matrix L = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.01, 1.5);    
    Matrix R = MatrixFactory::MakeRandomSparse(n, n, nnz, 0.01, 1.5);
    Matrix C = L + R;
    Matrix D = 2.5*C - R;
    
    std::cout << "\n\n";
    std::cout << "C = L + R produces\n" << C;
    std::cout << "D = 2.5*C - R produces\n" << D;


    return (0);
}
