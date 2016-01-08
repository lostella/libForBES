#include "IterativeSolver.h"

IterativeSolver::~IterativeSolver() {}

int IterativeSolver::run(){
	while(!stop()) {
		iterate();
	}
}

