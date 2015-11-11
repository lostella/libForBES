/* 
 * File:   OpDCT2.h
 * Author: Pantelis Sopasakis
 *
 * Created on September 15, 2015, 3:36 PM
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
 * The <a href="https://en.wikipedia.org/wiki/Discrete_cosine_transform#DCT-II">DCT-II
 * transform</a> is a linear square transformation \f$T:\mathbb{R}^n \to \mathbb{R}^n\f$
 * given by
 * 
 * \f[
 * (T(x))_k = \sum_{i=0}^{n} x_i \cos \left[ \frac{\pi}{n}\left(i+\frac{1}{2}\right)k \right]
 * \f]
 * 
 * for \f$k=0,\ldots, n-1\f$.
 * 
 * The discrete cosine transform, and especially this version of it - DCT-II - is popular in signal 
 * and image processing, especially for lossy compression.
 */
class OpDCT2 : public LinearOperator {
public:
    
    using LinearOperator::call;
    using LinearOperator::callAdjoint;
    
    OpDCT2();

    explicit OpDCT2(size_t n);

    virtual ~OpDCT2();

    virtual int call(Matrix& y, double alpha, Matrix& x, double gamma);
   
    virtual int callAdjoint(Matrix& y, double alpha, Matrix& x, double gamma);

    virtual std::pair<size_t, size_t> dimensionIn();

    virtual std::pair<size_t, size_t> dimensionOut();

    virtual bool isSelfAdjoint();

private:

    std::pair<size_t, size_t> m_dimension;

};

#endif	/* OPDCT2_H */

