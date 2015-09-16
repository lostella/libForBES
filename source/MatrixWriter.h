/* 
 * File:   MatrixWriter.h
 * Author: Pantelis Sopasakis
 *
 * Created on August 26, 2015, 7:55 PM
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

#ifndef MATRIXWRITER_H
#define	MATRIXWRITER_H

#include "Matrix.h"
#include "string.h"

#undef MATRIX_NROWS 
#undef MATRIX_NCOLS
#undef MATRIX_TYPE
#undef MATRIX_DATA

#define MATRIX_NROWS "nrows"
#define MATRIX_NCOLS "ncols"
#define MATRIX_TYPE  "type"
#define MATRIX_DATALENGTH "datalength"
#define MATRIX_DATA  "data"
#define MATRIX_NZ "nnz"
#define MATRIX_ENFORCE_DENSE_MODE "enforceDenseMode"

class MatrixWriter {
public:
    MatrixWriter(Matrix& m_matrix);
    virtual ~MatrixWriter();

    enum WriteFormat {
        PLAIN_TXT,
        JSON
    };

    void enforceDenseMode(bool t);

    void setWriteFormat(WriteFormat format);


    /**
     * Writes a matrix into a file.
     * @param fp
     */
    void write(FILE* fp);
private:

    Matrix& m_matrix;
    bool m_enforceDenseMode;
    WriteFormat m_format;

    void printJSON(FILE* fp);
    void printTXT(FILE* fp);

};

#endif	/* MATRIXWRITER_H */

