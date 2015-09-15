/* 
 * File:   OpDCT2.h
 * Author: chung
 *
 * Created on September 15, 2015, 3:36 PM
 */

#ifndef OPDCT2_H
#define	OPDCT2_H

#include "LinearOperator.h"
#include <math.h>

#define _USE_MATH_DEFINES


/**
 * \class OpDCT2
 * \brief The Discrete Cosine Transform Type II (DCT-II)
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on September 15, 2015, 3:36 PM
 * 
 * \ingroup LinOp
 * 
 * The DCT-II transformation is a linear square transformation \f$T:\mathbb{R}^n \to \mathbb{R}^n\f$
 * given by
 * 
 * \f[
 * (T(x))_k = \sum_{i=0}^{n} x_i \cos \left[ \frac{\pi}{n}(i+\frac{1}{2})k \right]
 * \f]
 */
class OpDCT2 : public LinearOperator {
public:
    OpDCT2();
    
    virtual ~OpDCT2();

    virtual Matrix call(Matrix& x);

    virtual Matrix callAdjoint(Matrix& x);

    virtual size_t dimensionIn();

    virtual size_t dimensionOut();

    virtual bool isSelfAdjoint();

private:
    

};

#endif	/* OPDCT2_H */

