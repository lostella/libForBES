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

    

    Matrix A = MatrixFactory::MakeRandomMatrix(20, 30, 0.0, 10.0, Matrix::MATRIX_DENSE);
    Matrix B = MatrixFactory::MakeRandomMatrix(10, 11, 0.0, 2.0, Matrix::MATRIX_DENSE);
    B.transpose();
    A.transpose();
        
    Matrix Asub = A.submatrixCopy(1, 3, 2, 5);  // 3 x 4
    Matrix Bsub = B.submatrixCopy(3, 6, 9, 10);  // 4 x 2
    Matrix exact = Asub*Bsub;
    

    std::cout << B;
    std::cout << Bsub;
    
    
    std::cout << "\n";
    Matrix result = A.multiplySubmatrix(B, 1,3,2,5, 3,6,9,10);
    std::cout << result-exact;

    return (0);
}


