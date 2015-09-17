/* 
 * File:   OpGradient2D.h
 * Author: chung
 *
 * Created on September 16, 2015, 6:20 PM
 */

#ifndef OPGRADIENT2D_H
#define	OPGRADIENT2D_H

#include "LinearOperator.h"


class OpGradient2D : public LinearOperator {
public:
    OpGradient2D();
    
    virtual ~OpGradient2D();

    virtual Matrix call(Matrix& x);

    virtual Matrix callAdjoint(Matrix& x);

    virtual size_t dimensionIn();

    virtual size_t dimensionOut();

    virtual bool isSelfAdjoint();


private:

};

#endif	/* OPGRADIENT2D_H */

