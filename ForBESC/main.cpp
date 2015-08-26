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
    
    std::cout << A;
    LDLFactorization * solver = new LDLFactorization(A);
    int status = solver->factorize();
    std::cout << "status = " << status << "\n";
    delete solver;

    //
    //    /* only the upper triangular part of A is required */
    //    int Ap [N + 1] = {0, 1, 2, 3, 4, 6, 7, 9, 11, 15, ANZ};
    //    int Ai[ANZ] =    {0, 1, 2, 3, 1, 4, 5, 4, 6, 4, 7, 0, 4, 7, 8, 1, 4, 6, 9};
    //    double Ax[ANZ] = {1.7, 1., 1.5, 1.1, .02, 2.6, 1.2, .16, 1.3, .09, 1.6,
    //        .13, .52, .11, 1.4, .01, .53, .56, 3.1};
    //    double b[N] = {.287, .22, .45, .44, 2.486, .72, 1.55, 1.424, 1.621, 3.759};
    //
    //    double Lx [LNZ], D [N], Y [N];
    //    int Li [LNZ], Lp [N + 1], Parent [N], Lnz [N], Flag [N], Pattern [N], d, i;
    //
    //    /* factorize A into LDL’ (P and Pinv not used) */
    //    ldl_symbolic(N, Ap, Ai, Lp, Parent, Lnz, Flag, NULL, NULL);
    //    printf("Non-zeros in L, excluding diagonal: %d\n", Lp [N]);
    //
    //    d = ldl_numeric(N, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Pattern,
    //            Flag, NULL, NULL);
    //    if (d == N) {
    //        /* solve Ax=b, overwriting b with the solution x */
    //        ldl_lsolve(N, b, Lp, Li, Lx);
    //        ldl_dsolve(N, b, D);
    //        ldl_ltsolve(N, b, Lp, Li, Lx);
    //        for (i = 0; i < N; i++) printf("x [%d] = %g\n", i, b [i]);
    //    } else {
    //        printf("ldl_numeric failed, D (%d,%d) is zero\n", d, d);
    //    }

    return (0);
}

