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

    OpLTI(Matrix& A, Matrix& B);

    virtual ~OpLTI();

    virtual Matrix call(Matrix& u);

    virtual Matrix callAdjoint(Matrix& x);

    virtual std::pair<size_t, size_t> dimensionIn();

    virtual std::pair<size_t, size_t> dimensionOut();

    virtual bool isSelfAdjoint();

private:

    Matrix & A;
    Matrix & B;

};

#endif	/* OPLTI_H */

