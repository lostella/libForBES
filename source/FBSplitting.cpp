#include "FBCache.h"
#include "FBSplitting.h"

int FBSplitting::prepare() {
  x = z;
  return 0;
}

int FBSplitting::iterate() {
  FBCache cache(prob, x, gamma);
  z = cache.forwardBackwardStep(gamma);
  return 0;
}

int FBSplitting::stop() {
  if (it >= maxIt) return 1;
  return 0;
}

FBSplitting::~FBSplitting() {}
