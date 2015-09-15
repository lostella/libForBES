/* 
 * File:   MatrixWriter.cpp
 * Author: Pantelis Sopasakis
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
    fprintf(fp, "  \"%s\":\"%s\",\n", MATRIX_TYPE, m_matrix.getTypeString().c_str());
    fprintf(fp, "  \"%s\":%zu,\n", MATRIX_NROWS, m_matrix.getNrows());
    fprintf(fp, "  \"%s\":%zu,\n", MATRIX_NCOLS, m_matrix.getNcols());
    if (m_matrix.getType() != Matrix::MATRIX_SPARSE) {
        fprintf(fp, "  \"%s\":%zu,\n", MATRIX_DATALENGTH, m_matrix.length());
    } else {
        fprintf(fp, "  \"%s\":%zu,\n", MATRIX_NZ, m_matrix.m_triplet->nnz);
    }

    fprintf(fp, "  \"%s\":%d,\n", MATRIX_ENFORCE_DENSE_MODE, m_enforceDenseMode);

    fprintf(fp, "  \"%s\":[", MATRIX_DATA);
    if (m_matrix.getType() == Matrix::MATRIX_DENSE || m_enforceDenseMode) {
        for (size_t j = 0; j < m_matrix.getNcols(); j++) {
            for (size_t i = 0; i < m_matrix.getNrows(); i++) {
                fprintf(fp, "%g", m_matrix.get(i, j));
                if (!(i == m_matrix.getNrows() - 1 && j == m_matrix.getNcols() - 1)) {
                    fprintf(fp, ", ");
                }
            }
        }
    } else if (m_matrix.getType() != Matrix::MATRIX_SPARSE) {
        for (size_t i = 0; i < m_matrix.length(); i++) {
            fprintf(fp, "%g", m_matrix.m_data[i]);
            if (i != m_matrix.length() - 1) {
                fprintf(fp, ", ");
            }
        }
    } else if (m_matrix.getType() == Matrix::MATRIX_SPARSE) {
        for (size_t i = 0; i < m_matrix.m_triplet->nnz; i++) {
            fprintf(fp, "[%d, %d, %g]", ((int*) m_matrix.m_triplet->i)[i], ((int*) m_matrix.m_triplet->j)[i], ((double*) m_matrix.m_triplet->x)[i]);
            if (i != m_matrix.m_triplet->nnz - 1) {
                fprintf(fp, ", ");
            }
        }
    }
    fprintf(fp, "]\n}");
}

void MatrixWriter::printTXT(FILE* fp) {
    fprintf(fp, "%d ", (int) (m_matrix.getType()));
    fprintf(fp, "%zu ", m_matrix.getNrows());
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

    if (m_format == JSON) {
        printJSON(fp);
    } else if (m_format = PLAIN_TXT) {
        printTXT(fp);
    }


}
