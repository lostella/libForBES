#include "FBCache.h"
#include "FBSplitting.h"

FBSplitting::FBSplitting(FBProblem & prob, Matrix & x0, double gamma) : m_cache(FBCache(prob, x0, gamma)) {
    m_prob = &prob;
    m_gamma = gamma;
    setMaxIt(1000);
    m_tol = 1e-6;
}

int FBSplitting::iterate() {
    Matrix * z = m_cache.get_forward_backward_step(m_gamma);
    m_cache.set_point(*z);
    return 0;
}

int FBSplitting::stop() {
    if (m_cache.get_norm_fpr() / m_gamma <= m_tol) return 1;
    return 0;
}

void FBSplitting::setTol(double tol) {
    m_tol = tol;
}

Matrix& FBSplitting::getSolution() {
    return *m_cache.get_forward_backward_step(m_gamma);
}

FBSplitting::~FBSplitting() {
	
}
