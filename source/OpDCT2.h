/* 
 * File:   OpDCT2.h
 * Author: Pantelis Sopasakis
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
 * (T(x))_k = \sum_{i=0}^{n} x_i \cos \left[ \frac{\pi}{n}\left(i+\frac{1}{2}\right)k \right]
 * \f]
 * 
 * for \f$k=0,\ldots, n-1\f$.
 * 
 * The DCT, and especially this version of it - DCT-II - is popular in signal 
 * and image processing, especially for lossy compression.
 */
class OpDCT2 : public LinearOperator {
public:
    OpDCT2();

    OpDCT2(size_t n);

    virtual ~OpDCT2();

    virtual Matrix call(Matrix& x);

    virtual Matrix callAdjoint(Matrix& x);

    virtual size_t dimensionIn();

    virtual size_t dimensionOut();

    virtual bool isSelfAdjoint();

private:

    size_t m_dimension;

};

#endif	/* OPDCT2_H */

