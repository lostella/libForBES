/* 
 * File:   SeparableSum.cpp
 * Author: chung
 * 
 * Created on October 30, 2015, 7:22 PM
 */

#include "SeparableSum.h"
#include <iostream>

SeparableSum::~SeparableSum() {
}

SeparableSum::SeparableSum(std::map<Function*, std::vector<size_t>*> fun_idx_map) :
Function(), m_fun_idx_map(fun_idx_map) {
}

int SeparableSum::call(Matrix& x, double& f) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP

    /* iterator for the map of functions and indices */
    std::map<Function*, std::vector<size_t> * >::iterator map_iterator;

    /* current function */
    Function * c_fun = NULL;

    /* current vector of indices*/
    std::vector<size_t> * c_idx = NULL;

    /* sub-matrix */
    Matrix *c_x = NULL;

    /* initialization of result f(x) */
    f = 0.0;

    /* temporary value of sub-invocations */
    double f_temp;

    /* status of sub-invocations */
    int status;

    /* Invoke all functions */
    for (map_iterator = m_fun_idx_map.begin(); map_iterator != m_fun_idx_map.end(); map_iterator++) {
        /* Pointer to function      */
        c_fun = map_iterator->first;

        /* Function indices         */
        c_idx = map_iterator->second;

        /* initialize new sub-matrix */
        c_x = new Matrix(c_idx->size(), 1);

        /* iterator for the indices of c_fun */
        std::vector<size_t>::iterator idx_iterator;
        size_t k = 0;

        /* construct the sub-matrix  */
        for (idx_iterator = c_idx->begin(); idx_iterator != c_idx->end(); idx_iterator++) {
            c_x -> set(k, 0, x.get(*idx_iterator, 0));
            k++;
        }

        /* invoke sub-function on c_x, return f_temp */
        status = c_fun -> call(*c_x, f_temp);
        //LCOV_EXCL_START
        if (ForBESUtils::STATUS_OK != status) {
            delete c_x;
            c_x = NULL;
            return status;
        }
        //LCOV_EXCL_STOP
        f += f_temp;
    }
    delete c_x;
    c_x = NULL;
    return ForBESUtils::STATUS_OK;
}

int SeparableSum::call(Matrix& x, double& f, Matrix& grad) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

int SeparableSum::callProx(const Matrix& x, double gamma, Matrix& prox) {
    //LCOV_EXCL_START
    if (!x.isColumnVector()) {
        throw std::invalid_argument("x must be a column-vector");
    }
    //LCOV_EXCL_STOP

    std::map<Function*, std::vector<size_t> * >::iterator map_iterator;
    Function * c_fun = NULL;
    std::vector<size_t> * c_idx = NULL;
    Matrix *c_x = NULL;
    Matrix *c_prox = NULL;
    int status;

    for (map_iterator = m_fun_idx_map.begin(); map_iterator != m_fun_idx_map.end(); map_iterator++) {
        c_fun = map_iterator->first;
        c_idx = map_iterator->second;
        c_x = new Matrix(c_idx->size(), 1);
        c_prox = new Matrix(c_idx->size(), 1);
        std::vector<size_t>::iterator idx_iterator;
        size_t k = 0;
        for (idx_iterator = c_idx->begin(); idx_iterator != c_idx->end(); idx_iterator++) {
            c_x -> set(k, 0, x.get(*idx_iterator, 0));
            k++;
        }
        status = c_fun -> callProx(*c_x, gamma, *c_prox);
        //LCOV_EXCL_START
        if (ForBESUtils::STATUS_OK != status) {
            delete c_x;
            delete c_prox;
            return status;
        }
        //LCOV_EXCL_STOP
        k = 0;
        for (idx_iterator = c_idx->begin(); idx_iterator != c_idx->end(); idx_iterator++) {
            prox.set(*idx_iterator, 0, c_prox->get(k, 0));
            k++;
        }
    }
    delete c_x;
    delete c_prox;

    return ForBESUtils::STATUS_OK;
}

int SeparableSum::callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    return ForBESUtils::STATUS_UNDEFINED_FUNCTION;
}

FunctionOntologicalClass SeparableSum::category() {
    FunctionOntologicalClass meta("SeparableSum");
    meta.set_defines_f(true);
    meta.set_defines_grad(true);
    meta.set_defines_prox(true);
    return meta;
}


