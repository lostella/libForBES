/*
 * File:   IndSOC.h
 * Author: Lorenzo Stella
 *
 * Created on September 21, 2015, 10:07 AM
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
 
#ifndef INDSOC_H
#define	INDSOC_H

#include "Matrix.h"
#include "Function.h"

/*
 * The following class implements the indicator function of second order cones.
 * If n is the dimension of the considered space, then this is
 *
 *		SOC(n) = {(x, t) in R^n : x in R^(n-1) and norm(x) <= t}
 *
 * The dimension n must be given to the (only) constructor.
 */

class IndSOC : public Function {
public:

    IndSOC(int n);

    virtual ~IndSOC();
    
    virtual int call(Matrix& x, double& f);
    
    virtual int category();
    
    virtual int callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    virtual int callProx(const Matrix& x, double gamma, Matrix& prox);

protected:

private:

    int n;

};

#endif	/* INDSOC_H */

