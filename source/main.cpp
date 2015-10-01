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

       const size_t n = 14;
    const size_t m = 11;
    Matrix A = MatrixFactory::MakeRandomMatrix(m, n, 2.0, 10.0, Matrix::MATRIX_DENSE);
    std::cout << A;
    Matrix subA = A.submatrixCopy(1, 2, 4, 8);
        std::cout << subA;

    return (0);
}


