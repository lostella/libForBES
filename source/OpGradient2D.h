/* 
 * File:   OpGradient2D.h
 * Author: chung
 *
 * Created on September 16, 2015, 6:20 PM
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

#ifndef OPGRADIENT2D_H
#define	OPGRADIENT2D_H

#include "LinearOperator.h"


class OpGradient2D : public LinearOperator {
public:
    
    using LinearOperator::call;
    using LinearOperator::callAdjoint;
    
    OpGradient2D();
    
    virtual ~OpGradient2D();

    virtual int call(Matrix& y, double alpha, Matrix& x, double gamma);

    virtual int callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma);

    virtual std::pair<size_t, size_t> dimensionIn();

    virtual std::pair<size_t, size_t> dimensionOut();

    virtual bool isSelfAdjoint();


private:

};

#endif	/* OPGRADIENT2D_H */

