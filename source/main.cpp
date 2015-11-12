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
#include "CGSolver.h"


using namespace std;

#define N 3

void f11();

void f11() {
    Matrix::destroy_handle();
    std::cout << "bye!";
}

int main(int argc, char** argv) {

    double z_data[10] = {1.050000000000000, 1.070000000000000,1.090000000000000,1.110000000000000,1.130000000000000,1.150000000000000,1.170000000000000,1.190000000000000,1.210000000000000,1.230000000000000};
    Matrix z(10,1,z_data);
    OpDCT2 T(10);
    
    Matrix r = T.callAdjoint(z);
    std::cout << r;




    return EXIT_SUCCESS;
}

