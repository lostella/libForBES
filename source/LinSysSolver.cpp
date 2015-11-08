/* 
 * File:   LinSysSolver.cpp
 * Author: chung
 * 
 * Created on November 7, 2015, 3:47 PM
 */

#include "LinSysSolver.h"

LinSysSolver::LinSysSolver(Matrix& matrix) : m_matrix(matrix) {
    m_matrix_nrows = matrix.getNrows();
    m_matrix_ncols = matrix.getNcols();
    m_matrix_type = matrix.getType();
}

LinSysSolver::~LinSysSolver() {
}

