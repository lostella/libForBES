#include "Solver.h"

Solver::~Solver() {}

int FBSolver::setTolerance(double tol) {
  this->tol = tol;
}

int FBSolver::setMaxIters(int maxIt) {
  this->maxIt = maxIt;
}
