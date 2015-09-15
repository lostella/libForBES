/* 
 * File:   OpTotalVariation.h
 * Author: chung
 *
 * Created on September 15, 2015, 2:50 PM
 */

#ifndef OPTOTALVARIATION_H
#define	OPTOTALVARIATION_H

#include "LinearOperator.h"


/**
 * \class OpTotalVariation
 * \brief Total variation operator
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on September 15, 2015, 2:50 PM
 * 
 * \ingroup LinOp
 */
class OpTotalVariation : public LinearOperator {
public:
    OpTotalVariation();
    OpTotalVariation(const OpTotalVariation& orig);
    virtual ~OpTotalVariation();
private:

};

#endif	/* OPTOTALVARIATION_H */

