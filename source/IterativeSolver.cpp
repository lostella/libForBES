#include "IterativeSolver.h"

IterativeSolver::IterativeSolver(int maxit) {
    m_it = 0;
    m_maxit = maxit;
}

int IterativeSolver::run() {
    int status = ForBESUtils::STATUS_OK;
    while (m_it < m_maxit && !stop() && !ForBESUtils::is_status_error(status)) {
        status = iterate();
        m_it++;
    }
    return status;
}

int IterativeSolver::getIt() {
    return m_it;
}

IterativeSolver::~IterativeSolver() {
}