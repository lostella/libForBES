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
    
    void write(FILE* fp);
private:
    
    Matrix& m_matrix;

};

#endif	/* MATRIXWRITER_H */

