/* 
 * File:   MatrixSolver.h
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

#ifndef MATRIXSOLVER_H
#define	MATRIXSOLVER_H

#include "LinSysSolver.h"

/**
 * \class MatrixSolver
 * \brief Abstraction layer for linear systems 
 * \version version 0.0
 * \date Created on November 11, 2015, 3:44 AM
 * \author Pantelis Sopasakis
 * \ingroup LinSysSolver-group
 * 
 */ 
class MatrixSolver : public LinSysSolver {
public:

    virtual ~MatrixSolver();
private:

protected:

    /**
     * Instantiate a new LinSysSolver object given the system matrix by reference.
     * This constructor is protected, it is, therefore, only to be invoked by
     * subclasses of %LinSysSolver.
     * @param matrix the system matrix \f$A\f$
     */
    explicit MatrixSolver(Matrix& matrix);

    /**
     * Underlying matrix.
     */
    Matrix * m_matrix;
    /**
     * Type of the matrix
     */
    Matrix::MatrixType m_matrix_type;
    /**
     * Rows of the matrix
     */
    size_t m_matrix_nrows;
    /**
     * Columns of the matrix
     */
    size_t m_matrix_ncols;

};

#endif	/* MATRIXSOLVER_H */

