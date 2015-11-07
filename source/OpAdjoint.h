/* 
 * File:   OpAdjoint.h
 * Author: Pantelis Sopasakis
 *
 * Created on September 15, 2015, 2:57 PM
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

#ifndef OPADJOINT_H
#define	OPADJOINT_H

#include "LinearOperator.h"

/**
 * \class OpAdjoint
 * \brief Adjoint of a given linear operator
 * \version 0.1
 * \author Pantelis Sopasakis
 * \date Created on September 15, 2015, 2:57 PM
 * 
 * \ingroup LinOp
 */
class OpAdjoint : public LinearOperator {
public:

    explicit OpAdjoint(LinearOperator& op);

    virtual ~OpAdjoint();

    virtual Matrix call(Matrix& x);

    virtual Matrix callAdjoint(Matrix& x);

    virtual std::pair<size_t, size_t> dimensionIn();

    virtual std::pair<size_t, size_t> dimensionOut();

    virtual bool isSelfAdjoint();

private:

    LinearOperator& m_originalOperator;

};

#endif	/* OPADJOINT_H */

