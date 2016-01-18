#include "IterativeSolver.h"

IterativeSolver::~IterativeSolver() {
}

int IterativeSolver::run() {
    int status = ForBESUtils::STATUS_OK;
    while (m_it < m_maxit && !stop() && !ForBESUtils::is_status_error(status)) {
        status = iterate();
        m_it++;
    }
    return status;
}

void IterativeSolver::setMaxIt(int maxit) {
    m_maxit = maxit;
}

int IterativeSolver::getIt() {
    return m_it;
}
