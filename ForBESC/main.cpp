/* 
 * File:   main.cpp
 * Author: chung
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

#include <set>


using namespace std;

int main(int argc, char** argv) {
    const size_t n = 10;
    Matrix x = MatrixFactory::MakeSparse(n, 1, 0, Matrix::SPARSE_UNSYMMETRIC);
    Matrix result;
    result = x * x;

    Matrix y = MatrixFactory::MakeSparse(n, 1, 1, Matrix::SPARSE_UNSYMMETRIC);
    y.set(0, 0, 2.0);
    Matrix z = y*y;
    std::cout << z;

    x = y;
    std::cout << x;

    z = x+x;
    return (0);
}

