/* 
 * File:   LinearOperator.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 24, 2015, 5:05 PM
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

#ifndef LINEAROPERATOR_H
#define	LINEAROPERATOR_H

#include "Matrix.h"

/**
 * \class LinearOperator
 * \brief A linear operator T(x)
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on July 24, 2015, 5:05 PM
 * 
 * An interface for an arbitrary linear operator.
 * This is an abstract class which cannot be instantiated. 
 */
class LinearOperator {
public:

    /**
     * 
     * @param x
     * @return 
     */
    virtual Matrix call(Matrix& x) = 0;
    virtual Matrix callAdjoint(Matrix& x) = 0;
    virtual bool isSelfAdjoint() = 0;
    virtual size_t dimensionIn() = 0;
    virtual size_t dimensionOut() = 0;

    virtual ~LinearOperator();

protected:
    LinearOperator();
    LinearOperator(const LinearOperator& orig);

};

#endif	/* LINEAROPERATOR_H */

