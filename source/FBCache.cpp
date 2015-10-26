#include "FBCache.h"

FBCache::FBCache(FBProblem &p, Matrix & x, double gamma) {
	this->f1 = p.getf1();
	this->L1 = p.getL1();
	this->d1 = p.getd1();
	this->f2 = p.getf2();
	this->L2 = p.getL2();
	this->d2 = p.getd2();
	this->lin = p.getlin();
	this->g = p.getg();
	this->x = x;
	this->gamma = gamma;
}

double FBCache::evalf() {
	if (flag_evalf == 1) {
		return 0;
	}

	if (f1) {
		if (L1) res1x = L1->call(x);
		else res1x = x;
		if (d1) res1x += *d1;
		f1->call(res1x, f1x, gradf1x);
	}

	if (f2) {
		if (L2) res2x = L2->call(x);
		else res2x = x;
		if (d2) res2x += *d2;
		f2->call(res2x, f2x);
	}

	if (lin) {
		lin->transpose();
		linx = ((*lin)*x)[0];
		lin->transpose();
	}

	fx = f1x + f2x + linx;
	return fx;
}

Matrix& FBCache::forwardStep(double gamma) {
	if (flag_evalf == 0) {
		this->evalf();
	}
	if (flag_gradstep == 1) {
		if (gamma != this->gamma) {
			this->gamma = gamma;
			y = x - gamma*gradfx;
		}
	}
	gradfx = 0;
	if (f1) {
		if (L1) gradfx += L1->callAdjoint(gradf1x);
		else gradfx += gradf1x;
	}
	if (f2) {
		gradf2x = f2->call(x, f2x, gradf2x);
		if (L2) gradfx += L2->callAdjoint(gradf2x);
		else gradfx += gradf2x;
	}
	if (lin) {
		gradfx += (*lin);
	}
	y = x - gamma*gradfx;
	flag_gradstep = 1;
	return y;
}

Matrix& FBCache::forwardBackwardStep(double gamma) {
	double gamma0 = this->gamma;
	if (flag_gradstep == 0 || gamma != gamma0) {
		this->forwardStep(gamma);
	}
	if (flag_proxgradstep == 0 || gamma != gamma0) {
		g->callProx(x, gamma, z, gz);
	}
	flag_proxgradstep = 1;
	return z;
}

double FBCache::evalFBE(double gamma) {

}

Matrix& FBCache::gradFBE(double gamma) {

}
