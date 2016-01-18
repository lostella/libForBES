#include "FBStoppingRelative.h"

#include <cmath>
#include <iostream>

FBStoppingRelative::FBStoppingRelative(double tol) : FBStopping(tol) {

}

int FBStoppingRelative::stop(FBCache & c) {
    double normx = 0;
    Matrix * x = c.get_point();
    for (int i = 0; i < x->length(); i++) {
        double acc = (*x)[i];
        normx += acc*acc;
    }
    normx = sqrt(normx);
    if (c.get_norm_fpr() <= m_tol * (1 + normx)) return 1;
    else return 0;
}

FBStoppingRelative::~FBStoppingRelative() {

}
