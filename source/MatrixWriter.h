/* 
 * File:   MatrixWriter.h
 * Author: Pantelis Sopasakis
 *
 * Created on August 26, 2015, 7:55 PM
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

