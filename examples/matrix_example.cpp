#include "ForBES.h"
#include <iostream>

int main(int argc, char** argv) {

    // define sizes
    const size_t n = 5, m = 7;

    // Create a dense matrix 5-by-7 with all entries set to 0.0:
    Matrix A(n, m, Matrix::MATRIX_DENSE);

    // Set A(1,0) to a value:
    A.set(1, 0, 10.56);

    // retrieve the value of A(1,0):
    double value = A.get(1, 0);

    std::cout << "value = " << value << std::endl;



    return (0);
}
