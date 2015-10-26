/*
 * File:   FBCache.h
 * Author: Lorenzo Stella
 *
 * Created on October 2, 2015
 *
 * ForBES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ForBES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ForBES. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FBCACHE_H
#define	FBCACHE_H

#include "Matrix.h"
#include "FBProblem.h"

class FBCache {
public:

  FBCache(FBProblem &p, Matrix & x, double gamma);

	virtual ~FBCache() {}

	double evalf();
  double evalFBE(double gamma);
  Matrix& gradFBE(double gamma);
  Matrix& forwardBackwardStep(double gamma);

private:

	int	flag_evalf, flag_gradstep, flag_proxgradstep;

	Function * f1, * f2, * g;
  LinearOperator * L1, * L2;
  Matrix * d1, * d2, * lin;

	Matrix x, y, z;
	Matrix res1x, gradf1x;
	Matrix res2x, gradf2x;
  Matrix gradfx;
  double f1x, f2x, linx, fx, gz, gamma;

  Matrix& forwardStep(double gamma);

protected:

};

#endif /* FBCACHE_H */
