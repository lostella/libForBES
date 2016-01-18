#include "FBCache.h"
#include "FBSplittingFast.h"
#include "FBStopping.h"

#include <iostream>

FBSplittingFast::FBSplittingFast(FBProblem & prob, Matrix & x0, double gamma) :
FBSplitting(prob, x0, gamma) {
    m_previous = NULL;
}

FBSplittingFast::FBSplittingFast(FBProblem & prob, Matrix & x0, double gamma, FBStopping & sc) :
FBSplitting(prob, x0, gamma, sc) {
    m_previous = NULL;
}

FBSplittingFast::FBSplittingFast(FBProblem & prob, Matrix & x0, double gamma, int maxit) :
FBSplitting(prob, x0, gamma, maxit) {
    m_previous = NULL;
}

FBSplittingFast::FBSplittingFast(FBProblem & prob, Matrix & x0, double gamma, FBStopping & sc, int maxit) :
FBSplitting(prob, x0, gamma, sc, maxit) {
    m_previous = NULL;
}

int FBSplittingFast::iterate() {
    // store current point temporarily
    Matrix * temp = new Matrix(*m_cache.get_point());
    // extrapolate if not first iterate
    if (m_previous != NULL) {
        // y = x + k/(k+2) (x - x')
        //   = (2k+2)/(k+2) x - k/(k+2) x'
        Matrix::add(*m_cache.get_point(), -(1.0 * m_it) / (m_it + 2), *m_previous, (2.0 * m_it + 2.0) / (m_it + 2));
        // tell FBCache that the point has changed
        m_cache.reset();
        // delete previously allocated 'm_previous'
        delete m_previous;
    }
    // store m_previous point
    m_previous = temp;
    // execute FBS iteration
    return FBSplitting::iterate();
}

int FBSplittingFast::stop() {
    return FBSplitting::stop();
}

Matrix& FBSplittingFast::getSolution() {
    return FBSplitting::getSolution();
}

FBSplittingFast::~FBSplittingFast() {
    if (m_previous != NULL) {
        delete m_previous;
        m_previous = NULL;
    }
}
