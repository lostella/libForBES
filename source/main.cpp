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

void createMatFile(const char * fname, Matrix F) {
    FILE *f;
    char buff[150];
    sprintf(buff, "/home/chung/Documents/MATLAB/libforbes/matrices/%s.mtx", fname);
    std::cout << buff << "\n";
    f = fopen(buff, "w");
    MatrixWriter writer(F);
    writer.enforceDenseMode(false);
    writer.setWriteFormat(MatrixWriter::JSON);
    writer.write(f);
    fclose(f);
}

int main(int argc, char** argv) {
    std::cout << "libforbes v0.2\n";
    Matrix D = MatrixFactory::MakeRandomMatrix(3, 7, 0.0, 1.0, Matrix::MATRIX_DENSE);
    Matrix L = MatrixFactory::MakeRandomMatrix(5, 5, 0.0, 1.0, Matrix::MATRIX_LOWERTR);
    Matrix H = MatrixFactory::MakeRandomMatrix(5, 5, 0.0, 1.0, Matrix::MATRIX_SYMMETRIC);
    Matrix X = MatrixFactory::MakeRandomMatrix(5, 5, 0.0, 1.0, Matrix::MATRIX_DIAGONAL);
    Matrix S = MatrixFactory::MakeRandomSparse(20, 20, 12, 0.0, 1.0);
    
    createMatFile("FDense", D);
    createMatFile("FLower", L);
    createMatFile("FSymmetric", H);
    createMatFile("FDiagonal", X);
    createMatFile("FSparse", S);
    
    
    
    std::cout << D.getTypeString() << std::endl;
    std::cout << H.getTypeString() << std::endl;
    std::cout << S.getTypeString() << std::endl;
    std::cout << L.getTypeString() << std::endl;
    std::cout << X.getTypeString() << std::endl;

    return (0);
}


