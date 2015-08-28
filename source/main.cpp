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

#include <set>


using namespace std;

#define N 10 /* A is 10-by-10 */
#define ANZ 19 /* # of nonzeros on diagonal and upper triangular part of A */
#define LNZ 13 /* # of nonzeros below the diagonal of L */

int main(int argc, char** argv) {

    Matrix A = MatrixFactory::MakeSparse(10, 10, 19, Matrix::SPARSE_SYMMETRIC_L);
    A.set(0, 0, 1.7);
    A.set(0, 8, 0.13);
    A.set(1, 1, 1.0);
    A.set(1, 9, 0.01);
    A.set(2, 2, 1.5);
    A.set(3, 3, 1.1);
    A.set(4, 4, 2.6);
    A.set(5, 5, 1.2);
    A.set(6, 6, 1.3);
    A.set(6, 6, 1.3);
    A.set(7, 7, 1.6);
    A.set(8, 8, 1.4);
    A.set(9, 9, 3.1);
    A.set(1, 4, 0.02);
    A.set(4, 6, 0.16);
    A.set(4, 7, 0.09);
    A.set(4, 8, 0.52);
    A.set(4, 9, 0.53);
    A.set(6, 9, 0.56);
    A.set(7, 8, 0.11);
    
    FILE *f = fopen("./matrices/F000005436.txt", "w");
    MatrixWriter writer(A);
    //writer.enforceDenseMode(true);
    writer.setWriteFormat(MatrixWriter::PLAIN_TXT);
    writer.write(f);
    fclose(f);
   

    return (0);
}


