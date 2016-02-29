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

#include <algorithm>

#include "SumOfNorm2.h"
#include "MatrixFactory.h"
#include "Norm2.h"

SumOfNorm2::SumOfNorm2(size_t k) : Norm(), m_partition_length(k) {
    m_mu = 1.0;
    m_norm2 = new Norm2();
}

SumOfNorm2::SumOfNorm2(double mu, size_t k) : Norm(), m_mu(mu), m_partition_length(k) {
    m_norm2 = new Norm2();
}

SumOfNorm2::~SumOfNorm2() {
    if (m_norm2 != NULL) {
        delete m_norm2;
        m_norm2 = NULL;
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
    int status = ForBESUtils::STATUS_OK;
    for (size_t j = 0; j < n_chunks && status == ForBESUtils::STATUS_OK; j++) {
        double f_temp;
        Matrix x_part = MatrixFactory::ShallowVector(x, m_partition_length, j * m_partition_length);
        status = m_norm2->call(x_part, f_temp);
        f += f_temp;
    }
    f *= m_mu;
    return status;
}

int SumOfNorm2::callProx(Matrix& v, double gamma, Matrix& prox) {
    const size_t n = v.getNrows();
    if (n % m_partition_length != 0) {
        throw std::invalid_argument("Given vector cannot be partitioned");
    }

    const size_t n_chunks = n / m_partition_length;
    Matrix v_part = MatrixFactory::ShallowVector();
    
    int status = ForBESUtils::STATUS_OK;
    for (size_t j = 0; j < n_chunks && status == ForBESUtils::STATUS_OK; j++) {
        double norm_temp;
        v_part = MatrixFactory::ShallowVector(v, m_partition_length, j * m_partition_length);
        Matrix prox_part = MatrixFactory::ShallowVector(prox, m_partition_length, j * m_partition_length);
        status = m_norm2->call(v_part, norm_temp);
        double factor = (1 - gamma / norm_temp);
        for (size_t i = 0; i < m_partition_length; i++) {
            prox_part[i] = factor * v_part[i];
        }
    }
    return status;
}

int SumOfNorm2::callProx(Matrix& v, double gamma, Matrix& prox, double& f_at_prox) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

FunctionOntologicalClass SumOfNorm2::category() {
    return FunctionOntologyRegistry::function();
}

int SumOfNorm2::dualNorm(Matrix& x, double& norm) {
    /* this implementation is very much similar to #call */
    norm = 0.0;
    const size_t n = x.getNrows();
    if (n % m_partition_length != 0) {
        throw std::invalid_argument("Given vector cannot be partitioned");
    }
    /* number of chunks */
    const size_t n_chunks = n / m_partition_length;    
    for (size_t j = 0; j < n_chunks; j++) {
        double f_temp;
        Matrix x_part = MatrixFactory::ShallowVector(x, m_partition_length, j * m_partition_length);
        m_norm2->call(x_part, f_temp);
        norm = std::max(norm, f_temp);
    }
    norm /= m_mu;
    return ForBESUtils::STATUS_OK;
}


