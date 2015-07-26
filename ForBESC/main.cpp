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
    size_t n = 10;
    Matrix Q = MatrixFactory::MakeRandomSparse(n, n, 20, 0.0, 1.0);
    Matrix Eye = MatrixFactory::MakeIdentity(n, 10.0);

    Q += Eye;

    Matrix Qt(Q);       // Qt= Q
    Qt.transpose();     // Qt = Qt'
    Q += Qt;            // Q = Q + Qt
   
    Function *F = new Quadratic(Q); // 0.5*x'Qx

    Matrix x = MatrixFactory::MakeRandomSparse(n, 1, 5, 0.0, 1.0);
    
    std::cout << x;
    
    

    double f; // F(x)
    f = Q.quad(x);
    int status = F->call(x, f);


    Matrix grad;
    double f_star;
    status = F->callConj(x, f_star, grad);

    std::cout << "f_star = " << f_star;

    std::cout << grad;

    return (0);
}

