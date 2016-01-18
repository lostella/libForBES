#include "FBStopping.h"

FBStopping::FBStopping(double tol) {
	m_tol = tol;
}

int FBStopping::stop(FBCache & c) {
	if (c.get_norm_fpr() <= m_tol) return 1;
	else return 0;
}

FBStopping::~FBStopping() {

}
