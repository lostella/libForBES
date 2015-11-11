/* 
 * File:   OpReverseVector.h
 * Author: Pantelis Sopasakis
 *
 * Created on September 15, 2015, 12:57 PM
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

#ifndef OPREVERSEVECTOR_H
#define	OPREVERSEVECTOR_H

#include "LinearOperator.h"
#include <algorithm>
#include <iterator>


/**
 * \class OpReverseVector
 * \brief A linear operator <code>T(x)</code> which reverses the order of elements in <code>x</code>
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on September 15, 2015, 12:57 PM
 * 
 * \ingroup LinOp
 */
class OpReverseVector : public LinearOperator {
public:
    
    using LinearOperator::call;
    using LinearOperator::callAdjoint;
    
    OpReverseVector();

    explicit OpReverseVector(size_t n);

    virtual ~OpReverseVector();

    virtual int call(Matrix& y, double alpha, Matrix& x, double gamma);

    virtual int callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma);

    virtual std::pair<size_t, size_t> dimensionIn();

    virtual std::pair<size_t, size_t> dimensionOut();

    virtual bool isSelfAdjoint();


private:

    /**
     * Input and output dimension of this operator.
     */
    size_t m_vectorDim;
};

#endif	/* OPREVERSEVECTOR_H */

