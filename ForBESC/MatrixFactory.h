/* 
 * File:   MatrixFactory.h
 * Author: chung
 *
 * Created on July 12, 2015, 7:50 PM
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

#ifndef MATRIXFACTORY_H
#define	MATRIXFACTORY_H

#include "Matrix.h"

class MatrixFactory {
    

public:
    MatrixFactory();
    MatrixFactory(const MatrixFactory& orig);
    virtual ~MatrixFactory();

    /**
     * Returns a random matrix as an instance of <code>Matrix</code>.
     * @param nrows Number of rows
     * @param ncols Number of columns
     * @param offset Random data follow a statistical distribution of the form 
     * <code>offset + scaling_factor * r</code>, where <code>r</code> is a random
     * variable in [0,1].     
     * @param scale The aforementioned scaling factor.
     * @param type matrix type.
     * @return A random matrix object.
     */

    static Matrix MakeRandomMatrix(int nrows, int ncols, float offset, float scale, Matrix::MatrixType type);
    static Matrix MakeIdentity(int n, float alpha);
    static Matrix MakeSparse(int nrows, int ncols, int max_nnz, Matrix::SparseMatrixType stype, cholmod_common *c);
    static Matrix MakeSparse(int nrows, int ncols, int max_nnz, Matrix::SparseMatrixType stype);
    static Matrix MakeSparseSymmetric(int nrows, int ncols, int max_nnz);

private:

};

#endif	/* MATRIXFACTORY_H */

