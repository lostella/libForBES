#include "FBProblem.h"

FBProblem::FBProblem(Function * f1, LinearOperator * L1, Matrix * d1, Function * f2, LinearOperator * L2, Matrix * d2, Matrix * lin, Function * g) {
    this->f1 = f1;
    this->L1 = L1;
    this->d1 = d1;
    this->f2 = f2;
    this->L2 = L2;
    this->d2 = d2;
    this->lin = lin;
    this->g = g;
}

