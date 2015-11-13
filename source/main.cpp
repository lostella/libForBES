/* 
 * File:   main.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on July 7, 2015, 7:47 PM
 * 
 * ForBES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *  
 * ForBES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with ForBES. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <map>

#include "ForBES.h"
#include "CGSolver.h"


using namespace std;

#define N 3

void f11();

void f11() {
    Matrix::destroy_handle();
    std::cout << "bye!";
}

int main(int argc, char** argv) {


    const size_t n = 11;
    Matrix T(n, n);
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            T.set(k, i, std::cos(static_cast<double> (k) * M_PI * (static_cast<double> (i) + 0.5) / static_cast<double> (n)));
        }
    }

    std::cout << T;
    Matrix x(n, 1);

    for (size_t j = 0; j < n; j++) {
        x[j] = j + 1;
    }

    std::cout << x;
    
    Matrix y_correct = T*x;

    LinearOperator * dct2 = new OpDCT2(n);
    Matrix y = dct2->call(x);
    
    Matrix E = y_correct - y;
    std::cout << E;
    

    delete dct2;



    return EXIT_SUCCESS;
}

