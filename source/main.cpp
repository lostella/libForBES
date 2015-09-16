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
#include "OpGradient.h"

#include <set>


using namespace std;

int main(int argc, char** argv) {

    size_t n = 10;

    LinearOperator * op = new OpGradient(n);

    Matrix x(n, 1);
    Matrix y(n - 1, 1);

    for (size_t i = 0; i < n; i++) {
        x.set(i, 0, i + 1);
    }
    
    for (size_t i = 0; i < n-1; i++) {
        y.set(i, 0, 3*i + 1);
    }

    Matrix Tx = op->call(x);
    Matrix Tstar_y = op->callAdjoint(y);

    Matrix inner1 = y*Tx;
    Matrix inner2 = x*Tstar_y;
    Matrix err = inner1 - inner2;

    std::cout << x;
    std::cout << y;
    
    std::cout << Tx;
    std::cout << Tstar_y;
    
    std::cout << err;
    std::cout << inner1;
    std::cout << inner2;



    delete op;

    return (0);
}


