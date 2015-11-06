/* 
 * File:   SumOfNorm2.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on November 6, 2015, 1:23 AM
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

#include "SumOfNorm2.h"
#include "MatrixFactory.h"

SumOfNorm2::SumOfNorm2(size_t k) : Function(), m_partition_index(k) {
    m_mu = 1.0;    
}

SumOfNorm2::SumOfNorm2(double mu, size_t k) : Function(), m_mu(mu), m_partition_index(k) {
}

SumOfNorm2::~SumOfNorm2() {
}

int SumOfNorm2::call(Matrix& x, double& f) {
    f = 0.0;
    Matrix x_part = MatrixFactory::ShallowVector();
    return 0;
}


