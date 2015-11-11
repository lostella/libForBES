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
    
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}


