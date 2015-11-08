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
#include "Norm2.h"

SumOfNorm2::SumOfNorm2(size_t k) : Function(), m_partition_length(k) {
    m_mu = 1.0;
    m_norm2 = new Norm2();
}

SumOfNorm2::SumOfNorm2(double mu, size_t k) : Function(), m_mu(mu), m_partition_length(k) {
    m_norm2 = new Norm2(mu);
}

SumOfNorm2::~SumOfNorm2() {
    if (m_norm2 != NULL) {
        delete m_norm2;
    }
}

int SumOfNorm2::call(Matrix& x, double& f) {
    f = 0.0;
    const size_t n = x.getNrows();
    if (n % m_partition_length != 0) {
        throw std::invalid_argument("Given vector cannot be partitioned");
    }
    /* number of chunks */
    const size_t n_chunks = n / m_partition_length;
    Matrix x_part = MatrixFactory::ShallowVector();
    double f_temp;
    for (size_t j = 0; j < n_chunks; j++) {
        x_part = MatrixFactory::ShallowVector(x, m_partition_length, j * m_partition_length);
        m_norm2->call(x, f_temp);
        f += f_temp;
    }
    return ForBESUtils::STATUS_OK;
}

int SumOfNorm2::callProx(Matrix& x, double gamma, Matrix& prox) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

int SumOfNorm2::callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

FunctionOntologicalClass SumOfNorm2::category() {
    return FunctionOntologyRegistry::function();
}



