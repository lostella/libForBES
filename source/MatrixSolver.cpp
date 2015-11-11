/* 
 * File:   MatrixSolver.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on November 11, 2015, 3:44 AM
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

#include "MatrixSolver.h"


MatrixSolver::MatrixSolver(Matrix& matrix) :
LinSysSolver(), m_matrix(&matrix) {
    m_matrix_nrows = matrix.getNrows();
    m_matrix_ncols = matrix.getNcols();
    m_matrix_type = matrix.getType();
}

MatrixSolver::~MatrixSolver() {
}

