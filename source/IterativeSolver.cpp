#include "IterativeSolver.h"

IterativeSolver::~IterativeSolver() {}

int IterativeSolver::run(){
	while (m_it < m_maxit && !stop()) {
		iterate();
		m_it++;
	}
}

void IterativeSolver::setMaxIt(int maxit) {
    m_maxit = maxit;
}

int IterativeSolver::getIt() {
    return m_it;
}
