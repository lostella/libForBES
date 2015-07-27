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

#include <set>


using namespace std;

int main(int argc, char** argv) {
    size_t n = 10;
   
    Matrix A = MatrixFactory::MakeRandomMatrix(n,n,0.0, 1.0, Matrix::MATRIX_LOWERTR);
    Matrix B = A;
    B.transpose();
    B *= -10.0;    

    std::cout << A;
    std::cout << B;
    
    return (0);
}

