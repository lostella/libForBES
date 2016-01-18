#include "IterativeSolver.h"

IterativeSolver::IterativeSolver(int maxit) {
	m_it = 0;
	m_maxit = maxit;
}

int IterativeSolver::run(){
	while (m_it < m_maxit && !stop()) {
		iterate();
		m_it++;
	}
}

int IterativeSolver::getIt() {
    return m_it;
}

IterativeSolver::~IterativeSolver() {}
