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
    cholmod_sparse *A, *B, *C;


    double alpha[1] = {1};
    double beta[1] = {1};
    cholmod_common c;
    cholmod_common c2;
    cholmod_start(&c); /* start CHOLMOD */
    cholmod_start(&c2); /* start CHOLMOD */

    FILE *fp;
    FILE *fp2;
    fp = fopen("matrices/sparse2.mx", "r");
    fp2 = fopen("matrices/sparse3.mx", "r");

    A = cholmod_read_sparse(fp, &c); /* read in a matrix */
    B = cholmod_read_sparse(fp2, &c2); /* read in a matrix */
    
    C = cholmod_add(A, B, alpha, beta, true, true, &c2);
    
    std::cout << C->nzmax << std::endl;
    std::cout << ((double*)C->x)[0];

    return (0);
}

