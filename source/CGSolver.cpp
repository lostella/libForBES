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
#include <iostream>
#include <lapacke.h>

CGSolver::CGSolver(LinearOperator& linop) : LinOpSolver(linop) {
    m_precond = NULL;
}

CGSolver::CGSolver(LinearOperator& linop, LinearOperator& preconditioner) : LinOpSolver(linop), m_precond(&preconditioner) {

}

CGSolver::~CGSolver() {
}

int CGSolver::solve(Matrix& b, Matrix& solution) const {
    Matrix r = b;
    Matrix z = m_precond->call(r);
    Matrix p = z;
    size_t k = 0;   
    bool keepgoing = true;
    Matrix Ap;
    while (keepgoing) {
        Ap = m_linop->call(p);                      // Ap = A * p
        double a_denom = (p * Ap).get(0, 0);        //
        double a_numer = (r * z).get(0, 0);         //
        double alpha = a_numer / a_denom;           // alpha = (r,z)/(p, Ap);
        Matrix::add(solution, alpha, p, 1.0);       // x = x + alpha * p
        Matrix r_new = r;                           // r_new = r
        double err = LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', r_new.getNrows(), 1, r_new.getData(), r_new.getNrows());
        std::cout << "*" << err << std::endl; 
        if (err < 1e-4){
            break;
        }
        Matrix::add(r_new, -alpha, Ap, 1.0);        // r_new = r - alpha A p;
        Matrix z_new = m_precond->call(r_new);      // z_new = P(r_new)
        double zr = (r_new * z_new).get(0, 0);      
        double beta = zr / a_numer;
        Matrix::add(p, 1.0, z_new, beta);
        z = z_new;
        r = r_new;
        k++;
        if (k > 1000) {
            keepgoing = false;
        }
    }
    std::cout << "solved in " << k << " iterations\n";
    return ForBESUtils::STATUS_OK;
}


