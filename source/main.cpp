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
#include "OpDCT3.h"

#include <set>


using namespace std;

int main(int argc, char** argv) {

    size_t n = 8;

    LinearOperator * op = new OpDCT3(n);

    Matrix x(n, 1);
    Matrix y(n, 1);
    

    for (size_t i = 0; i < n; i++) {
        x.set(i, 0, 0.9*i + 1.1);
        y.set(i, 0, 3*i + 0.5);
    }
        
    Matrix Tx = op->call(x);
    Matrix Tstar_y = op->callAdjoint(y);

    Matrix err = y*Tx - x*Tstar_y;
    
    std::cout << err;
    

    delete op;

    return (0);
}


