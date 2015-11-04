/* 
 * File:   OpLinearCombination.h
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

#ifndef OPLINEARCOMBINATION_H
#define	OPLINEARCOMBINATION_H

#include "LinearOperator.h"



/**
 * \class OpLinearCombination
 * \brief Linear combination of two linear operators <code>T(x) = a*A(x) + b*B(x)</code>
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on September 14, 2015, 9:25 PM
 * 
 * \ingroup LinOp
 */
class OpLinearCombination : public LinearOperator {
    
public:
    
    OpLinearCombination(LinearOperator& A, LinearOperator& B, double a, double b);

    virtual ~OpLinearCombination();
private:
    LinearOperator& m_A;
    LinearOperator& m_B;
    double m_a;
    double m_b;
};

#endif	/* OPLINEARCOMBINATION_H */

