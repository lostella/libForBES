/* 
 * File:   OpLinearCombination.h
 * Author: chung
 *
 * Created on September 14, 2015, 9:25 PM
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
    OpLinearCombination();
    OpLinearCombination(const OpLinearCombination& orig);
    virtual ~OpLinearCombination();
private:

};

#endif	/* OPLINEARCOMBINATION_H */

