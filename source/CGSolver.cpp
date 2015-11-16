/* 
 * File:   CGSolver.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on November 11, 2015, 12:41 AM
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

#include "CGSolver.h"
#include <lapacke.h>
#include <math.h>

CGSolver::CGSolver(LinearOperator& linop) : LinOpSolver(linop) {
    init();
    m_precond = NULL;
}

CGSolver::CGSolver(LinearOperator& linop, LinearOperator& preconditioner) : LinOpSolver(linop), m_precond(&preconditioner) {
    init();
}

CGSolver::~CGSolver() {
}

void CGSolver::init() {
    m_tolerance = 1e-4;
    m_max_iterations = 500;
    m_err = NAN;
    m_iterations_count = 0;
}

int CGSolver::solve(Matrix& b, Matrix& solution) {
    m_iterations_count = 0;
    m_err = NAN;
    Matrix r(b);
    Matrix z = m_precond->call(r);
    Matrix p = MatrixFactory::ShallowVector(z, 0);
    bool keepgoing = true;
    Matrix Ap;
    while (keepgoing) {
        Ap = m_linop->call(p);                      // Ap = A * p
        double a_denom = (p * Ap).get(0, 0);        //
        double a_numer = (r * z).get(0, 0);         //
        double alpha = a_numer / a_denom;           // alpha = (r,z)/(p, Ap);
        Matrix::add(solution, alpha, p, 1.0);       // x = x + alpha * p
        Matrix r_new = r;                           // r_new = r
        
        m_err = LAPACKE_dlange(
                LAPACK_COL_MAJOR,
                'I',
                r_new.getNrows(),
                1,
                r_new.getData(),
                r_new.getNrows());                  // compute ||r||_{inf}
        
        if (m_err < m_tolerance) {
            break;                                  // stop if tolerance reached
        }
        
        Matrix::add(r_new, -alpha, Ap, 1.0);        // r_new = r - alpha A p;
        Matrix z_new = m_precond->call(r_new);      // z_new = P(r_new)
        double zr = (r_new * z_new).get(0, 0);      
        double beta = zr / a_numer;                 // beta = (r_new, z_nwq)/(r, z)
        Matrix::add(p, 1.0, z_new, beta);           // p = beta p + z_new
        z = z_new;                                  // z = z_new
        r = r_new;                                  // r = r_new
        m_iterations_count++;                       // k = k + 1
        if (m_iterations_count > m_max_iterations) {
            keepgoing = false;
        }
    }
    return ForBESUtils::STATUS_OK;
}

double CGSolver::last_error() const {
    return m_err;
}

size_t CGSolver::last_num_iter() const {
    return m_iterations_count;
}


