#include "FBCache.h"
#include "FBSplitting.h"

FBSplitting::FBSplitting(FBProblem & prob, Matrix & x0, double gamma) : m_cache(FBCache(prob, x0, gamma)) {
    m_prob = &prob;
    m_gamma = gamma;
    m_it = 0;
    m_maxit = 1000;
}

int FBSplitting::iterate() {
    Matrix * z = m_cache.get_forward_backward_step(m_gamma);
    m_cache.set_point(*z);
    m_it++;
    return 0;
}

int FBSplitting::stop() {
    if (m_it >= m_maxit) return 1;
    if (m_cache.get_norm_fpr() / m_gamma <= m_tol) return 1;
    return 0;
}

int FBSplitting::setMaxIt(int maxit) {
    m_maxit = maxit;
}

int FBSplitting::setTol(double tol) {
    m_tol = tol;
}

int FBSplitting::getIt() {
    return m_it;
}

Matrix& FBSplitting::getSolution() {
    return *m_cache.get_forward_backward_step(m_gamma);
}

FBSplitting::~FBSplitting() {
}
