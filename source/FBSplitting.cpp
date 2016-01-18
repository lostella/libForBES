#include "FBCache.h"
#include "FBSplitting.h"
#include "FBStopping.h"

#include <iostream>

#define DEFAULT_MAXIT 1000
#define DEFAULT_TOL 1e-6

FBSplitting::FBSplitting(FBProblem & prob, Matrix & x0, double gamma) :
IterativeSolver(DEFAULT_MAXIT), m_cache(FBCache(prob, x0, gamma)) {
    m_prob = &prob;
    m_gamma = gamma;
    m_sc = new FBStopping(DEFAULT_TOL);
    delete_sc = true;
}

FBSplitting::FBSplitting(FBProblem & prob, Matrix & x0, double gamma, FBStopping & sc) :
IterativeSolver(DEFAULT_MAXIT), m_cache(FBCache(prob, x0, gamma)) {
    m_prob = &prob;
    m_gamma = gamma;
    m_sc = &sc;
    delete_sc = false;
}

FBSplitting::FBSplitting(FBProblem & prob, Matrix & x0, double gamma, int maxit) :
IterativeSolver(maxit), m_cache(FBCache(prob, x0, gamma)) {
    m_prob = &prob;
    m_gamma = gamma;
    m_sc = new FBStopping(DEFAULT_TOL);
    delete_sc = true;
}

FBSplitting::FBSplitting(FBProblem & prob, Matrix & x0, double gamma, FBStopping & sc, int maxit) :
IterativeSolver(maxit), m_cache(FBCache(prob, x0, gamma)) {
    m_prob = &prob;
    m_gamma = gamma;
    m_sc = &sc;
    delete_sc = false;
}

int FBSplitting::iterate() {
    m_cache.set_point(*m_cache.get_forward_backward_step(m_gamma));
    return 0;
}

int FBSplitting::stop() {
    return m_sc->stop(m_cache);
}

Matrix& FBSplitting::getSolution() {
    return *m_cache.get_forward_backward_step(m_gamma);
}

FBSplitting::~FBSplitting() {
    if (delete_sc && m_sc != NULL) {
        delete m_sc;
        m_sc = NULL;
    }
}
