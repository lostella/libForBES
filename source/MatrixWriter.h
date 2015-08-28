/* 
 * File:   MatrixWriter.h
 * Author: chung
 *
 * Created on August 26, 2015, 7:55 PM
 */

#ifndef MATRIXWRITER_H
#define	MATRIXWRITER_H

#include "Matrix.h"

class MatrixWriter {
public:
    MatrixWriter(Matrix& m_matrix);
    virtual ~MatrixWriter();
    
    enum WriteFormat{
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

