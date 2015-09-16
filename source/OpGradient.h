/* 
 * File:   OpGradient.h
 * Author: chung
 *
 * Created on September 16, 2015, 1:35 AM
 */

#ifndef OPGRADIENT_H
#define	OPGRADIENT_H

#include "LinearOperator.h"

class OpGradient : public LinearOperator {
public:
    OpGradient();

    OpGradient(size_t n);

    virtual ~OpGradient();

    virtual Matrix call(Matrix& x);

    virtual Matrix callAdjoint(Matrix& x);

    virtual size_t dimensionIn();

    virtual size_t dimensionOut();

    virtual bool isSelfAdjoint();


private:

    size_t m_dimension;

};

#endif	/* OPGRADIENT_H */

