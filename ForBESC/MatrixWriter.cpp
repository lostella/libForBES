/* 
 * File:   MatrixWriter.cpp
 * Author: chung
 * 
 * Created on August 26, 2015, 7:55 PM
 */

#include "MatrixWriter.h"

MatrixWriter::MatrixWriter(Matrix& m_matrix) : m_matrix(m_matrix) {
    m_enforceDenseMode = false;
    m_format = PLAIN_TXT;
}

MatrixWriter::~MatrixWriter() {
}

void MatrixWriter::enforceDenseMode(bool t) {
    m_enforceDenseMode = t;
}

void MatrixWriter::setWriteFormat(WriteFormat format) {
    m_format = format;
}

void MatrixWriter::printJSON(FILE* fp) {
    fprintf(fp, "{\n");
    fprintf(fp, "\t\"%s\" = %d,\n", "type", (int) (m_matrix.getType()));
    fprintf(fp, "\t\"%s\" = %zu,\n", "rows", m_matrix.getNrows());
    fprintf(fp, "\t\"%s\" = %zu,\n", "cols", m_matrix.getNcols());
    fprintf(fp, "\t\"data\" = [");

    if (m_matrix.getType() == Matrix::MATRIX_DENSE || m_enforceDenseMode) {
        for (size_t j = 0; j < m_matrix.getNcols(); j++) {
            for (size_t i = 0; i < m_matrix.getNrows(); i++) {
                fprintf(fp, "%g", m_matrix.get(i, j));
                if (i != m_matrix.getNrows() - 1 && j != m_matrix.getNcols() - 1) {
                    fprintf(fp, ", ");
                }
            }
            fprintf(fp, ", ");
        }
    }
    fprintf(fp, "]\n}");
}

void MatrixWriter::printTXT(FILE* fp) {
    fprintf(fp, "%d\n", (int) (m_matrix.getType()));
    fprintf(fp, "%zu\n", m_matrix.getNrows());
    fprintf(fp, "%zu\n", m_matrix.getNcols());

    if (m_matrix.getType() == Matrix::MATRIX_DENSE || m_enforceDenseMode) {
        for (size_t j = 0; j < m_matrix.getNcols(); j++) {
            for (size_t i = 0; i < m_matrix.getNrows(); i++) {
                fprintf(fp, "%g\n", m_matrix.get(i, j));
            }
        }
    }
}

void MatrixWriter::write(FILE* fp) {
    if (fp == NULL) {
        throw std::invalid_argument("File is NULL");
    }
    if (m_enforceDenseMode || m_matrix.getType() != Matrix::MATRIX_SPARSE) {
        if (m_format == JSON) {
            printJSON(fp);
        } else if (m_format == PLAIN_TXT) {
            printTXT(fp);
        }
    } else if (m_matrix.getType() == Matrix::MATRIX_SPARSE) {
        throw std::logic_error("Sparse serialization - not implemented yet!");
    }

}
