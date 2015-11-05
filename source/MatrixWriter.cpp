/* 
 * File:   MatrixWriter.cpp
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

#include "MatrixWriter.h"

MatrixWriter::MatrixWriter(Matrix& matrix) : m_matrix(matrix) {
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
    fprintf(fp, "  \"%s\":"FMT_SIZE_T",\n", MATRIX_NROWS, m_matrix.getNrows());
    fprintf(fp, "  \"%s\":"FMT_SIZE_T",\n", MATRIX_NCOLS, m_matrix.getNcols());
    if (m_matrix.getType() != Matrix::MATRIX_SPARSE) {
        fprintf(fp, "  \"%s\":"FMT_SIZE_T",\n", MATRIX_DATALENGTH, m_matrix.length());
    } else {
        fprintf(fp, "  \"%s\":"FMT_SIZE_T",\n", MATRIX_NZ, m_matrix.m_triplet->nnz);
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
            fprintf(fp, "[%d, %d, %g]",
                    (static_cast<int*> (m_matrix.m_triplet->i))[i],
                    (static_cast<int*> (m_matrix.m_triplet->j))[i],
                    (static_cast<double*> (m_matrix.m_triplet->x))[i]);
            if (i != m_matrix.m_triplet->nnz - 1) {
                fprintf(fp, ", ");
            }
        }
    }
    fprintf(fp, "]\n}");
}

void MatrixWriter::printTXT(FILE* fp) {
    fprintf(fp, "%d ", static_cast<int> (m_matrix.getType()));
    fprintf(fp, FMT_SIZE_T" ", m_matrix.getNrows());
    fprintf(fp, FMT_SIZE_T"\n", m_matrix.getNcols());

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
    } else if (m_format == PLAIN_TXT) {
        printTXT(fp);
    }


}
