/* 
 * File:   SumOfNorm2.h
 * Author: Pantelis Sopasakis
 *
 * Created on November 6, 2015, 1:23 AM
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

#ifndef SUMOFNORM2_H
#define	SUMOFNORM2_H

#include "Function.h"



/**
 * \class SumOfNorm2
 * \brief Sum of 2-norms
 * \version version 0.1
 * \ingroup Functions
 * \date Created on November 6, 2015, 1:23 AM
 * \author Pantelis Sopasakis
 * 
 */ 
class SumOfNorm2 : public Function {
public:
    using Function::call;
    
    SumOfNorm2(size_t k);
    
    SumOfNorm2(double mu, size_t k);

    virtual ~SumOfNorm2();
    
    virtual int call(Matrix& x, double& f);


private:
    double m_mu;
    size_t m_partition_index;

};

#endif	/* SUMOFNORM2_H */

