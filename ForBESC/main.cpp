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

    FILE *fp;
    fp = fopen("sparse1.mx", "r");

    if (fp == NULL) {
        fprintf(stderr, "Can't open input file in.list!\n");
        exit(1);
    }

    Matrix A = MatrixFactory::ReadSparse(fp);
    cout << A;

    return (0);
}

