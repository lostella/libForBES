/* 
 * File:   OpLinearCombination.h
 * Author: Pantelis Sopasakis
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
    
    OpLinearCombination(LinearOperator& A, LinearOperator& B, double a, double b);

    virtual ~OpLinearCombination();
private:
    LinearOperator& A;
    LinearOperator& B;
    double a;
    double b;
};

#endif	/* OPLINEARCOMBINATION_H */

