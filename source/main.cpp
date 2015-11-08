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
#include <stdio.h>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <map>

#include "ForBES.h"
#include "FunctionOntologyRegistry.h"
#include "LogLogisticLoss.h"
#include "Norm1.h"
#include "S_LDLFactorization.h"
#include "ConjugateFunction.h"


using namespace std;

#define N 3

int main(int argc, char** argv) {

    Matrix b = MatrixFactory::MakeRandomMatrix(5, 1, -1.0, 2.0);
    for (size_t i = 0; i < 5; i++) {
        b.set(i, 0, 0.9 * i + 1.0);
    }


    Function * f = new HingeLoss(b);
    Function * f_conj = new ConjugateFunction(*f);

    Matrix x = MatrixFactory::MakeRandomMatrix(5, 1, -1.0, 2.0);
    for (size_t i = 0; i < 5; i++) {
        x.set(i, 0, 0.1 * i + 0.08);
    }

    std::cout << x;

    Matrix prox(5, 1);
    f->callProx(x, 0.5, prox);

    std::cout << prox;
    prox = Matrix(5, 1);

    f_conj ->callProx(x, 0.5, prox);
    std::cout << prox;

    delete f;
    delete f_conj;

    return (0);
}

