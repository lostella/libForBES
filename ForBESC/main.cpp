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

using namespace std;

void add_A_entry(cholmod_triplet* A, int r, int c, double x) {
    ((int*) A->i)[A->nnz] = r;
    ((int*) A->j)[A->nnz] = c;
    ((double*) A->x)[A->nnz] = x;
    (A->nnz)++;
}

int main(int argc, char** argv) {

    cholmod_triplet *A_triplets;
    cholmod_sparse *A;
    cholmod_dense *x, *b, *r;
    cholmod_factor *L;
    double one [2] = {1, 0}, m1 [2] = {-1, 0}; /* basic scalars */
    cholmod_common c;
    cholmod_start(&c); /* start CHOLMOD */

    A_triplets = cholmod_allocate_triplet(3, 3, 4, 1, CHOLMOD_REAL, &c);
    add_A_entry(A_triplets, 0, 0, 4.0);
    add_A_entry(A_triplets, 0, 1, 1.0);
    //add_A_entry(A_triplets, 1, 0, 1.0);
    add_A_entry(A_triplets, 1, 1, 5.0);
    add_A_entry(A_triplets, 2, 2, 10.0);

    A = cholmod_triplet_to_sparse(A_triplets, 4, &c);


    b = cholmod_ones(A->nrow, 1, A->xtype, &c); /* b = ones(n,1) */
    L = cholmod_analyze(A, &c); /* analyze */
    cholmod_factorize(A, L, &c); /* factorize */
    x = cholmod_solve(CHOLMOD_A, L, b, &c); /* solve Ax=b */
    r = cholmod_copy_dense(b, &c); /* r = b */
    cholmod_sdmult(A, 0, m1, one, x, r, &c); /* r = r-Ax */
    printf("norm(b-Ax) %8.1e\n",
            cholmod_norm_dense(r, 0, &c)); /* print norm(r) */

    cholmod_sparse *Ls;
    cholmod_triplet *Lt;
    Ls = cholmod_factor_to_sparse(L, &c);
    Lt = cholmod_sparse_to_triplet(Ls, &c);

    std::cout << ((double*) Lt->x)[0] << "\n";
    std::cout << ((double*) Lt->x)[1] << "\n";
    std::cout << ((double*) Lt->x)[2] << "\n";
    std::cout << ((double*) Lt->x)[3] << "\n";

    cholmod_free_factor(&L, &c); /* free matrices */
    cholmod_free_sparse(&Ls, &c); /* free matrices */
    cholmod_free_triplet(&Lt, &c); /* free matrices */
    cholmod_free_sparse(&A, &c);
    cholmod_free_dense(&r, &c);
    cholmod_free_dense(&x, &c);
    cholmod_free_dense(&b, &c);

    cholmod_finish(&c); /* finish CHOLMOD */
    return (0);
}

