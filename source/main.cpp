/* 
 * File:   main.cpp
 * Author: Pantelis Sopasakis
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

#include "ldl.h"
#include "MatrixWriter.h"
#include "MatrixOperator.h"
#include "OpComposition.h"
#include "OpReverseVector.h"

#include <set>


using namespace std;

int main(int argc, char** argv) {

    LinearOperator * rev = new OpReverseVector();
    
    
    Matrix x = MatrixFactory::MakeRandomMatrix(5, 1, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix xrev = rev->call(x);
    
    std::cout << x;
    std::cout << xrev;
    
    return (0);
}


