/* 
 * File:   OpSum.h
 * Author: chung
 *
 * Created on September 14, 2015, 9:25 PM
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