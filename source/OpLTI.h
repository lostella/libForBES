/* 
 * File:   OpLTI.h
 * Author: Pantelis Sopasakis
 *
 * Created on September 30, 2015, 6:20 PM
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

#ifndef OPLTI_H
#define	OPLTI_H

#include "LinearOperator.h"

/**
 * \class OpLTI
 * \brief %OpLTI simulates a LTI system with zero initial condition
 * \version 0.0
 * \author Pantelis Sopasakis
 * \date Created on September 30, 2015, 6:20 PM
 * 
 * \ingroup LinOp
 */
class OpLTI : public LinearOperator {
public:
    
    using LinearOperator::call;
    using LinearOperator::callAdjoint;

    OpLTI(Matrix& A, Matrix& B, size_t N_horizon);

    virtual ~OpLTI();

    virtual int call(Matrix& y, double alpha, Matrix& x, double gamma);

    virtual int callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma);
    
    virtual std::pair<size_t, size_t> dimensionIn();

    virtual std::pair<size_t, size_t> dimensionOut();

    virtual bool isSelfAdjoint();

private:

    Matrix & m_A;
    Matrix & m_B;
    size_t m_N;

};

#endif	/* OPLTI_H */

