/* 
 * File:   OpSum.h
 * Author: Pantelis Sopasakis
 *
 * Created on September 14, 2015, 9:25 PM
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

#ifndef OPSUM_H
#define	OPSUM_H

#include "LinearOperator.h"


/**
 * \class OpSum
 * \brief The sum of two linear operators <code>T(x) = A(x) + B(x)</code>
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on September 14, 2015, 9:25 PM
 * 
 * \ingroup LinOp
 */
class OpSum : public LinearOperator {
public:
    OpSum();
    OpSum(const OpSum& orig);
    virtual ~OpSum();
private:

};

#endif	/* OPSUM_H */