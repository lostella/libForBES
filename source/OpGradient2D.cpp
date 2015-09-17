/* 
 * File:   OpGradient2D.cpp
 * Author: chung
 * 
 * Created on September 16, 2015, 6:20 PM
 */

#include "OpGradient2D.h"

OpGradient2D::OpGradient2D() {
}

OpGradient2D::~OpGradient2D() {
}

Matrix OpGradient2D::call(Matrix& x) {
    throw std::logic_error("NIY");
}

Matrix OpGradient2D::callAdjoint(Matrix& x) {
    throw std::logic_error("NIY");
}

size_t OpGradient2D::dimensionIn() {
    throw std::logic_error("NIY");
}

size_t OpGradient2D::dimensionOut() {
    throw std::logic_error("NIY");
}

bool OpGradient2D::isSelfAdjoint() {
    throw std::logic_error("NIY");
}

