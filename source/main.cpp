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
#include "OpDCT2.h"

#include <set>


using namespace std;

int main(int argc, char** argv) {

    size_t n = 10;

    LinearOperator * op = new OpDCT2();

    Matrix x(n, 1);
    

    for (size_t i = 0; i < n; i++) {
        x.set(i, 0, i + 1.0);
    }
        
    Matrix Tx = op->call(x);
    
    
    std::cout << x;
    std::cout << Tx;
    

    delete op;

    return (0);
}


