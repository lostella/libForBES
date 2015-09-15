/* 
 * File:   OpComposition.h
 * Author: chung
 *
 * Created on September 14, 2015, 9:26 PM
 */

#ifndef OPCOMPOSITION_H
#define	OPCOMPOSITION_H

#include "LinearOperator.h"

/**
 * \class OpComposition
 * \brief The composition of two linear operators <code>T(x) = A(B(x))</code>
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on September 14, 2015, 9:26 PM
 * 
 * \ingroup LinOp
 */
class OpComposition : public LinearOperator {
public:        

    OpComposition(LinearOperator& A, LinearOperator& B) : LinearOperator(), A(A), B(B) {
        // check dimensions
        if (A.dimensionIn() != B.dimensionOut()) {
            throw std::invalid_argument("A and B have incompatible dimensions; AoB is not well defined.");
        }
    }

    virtual ~OpComposition();

    virtual Matrix call(Matrix& x);

    virtual size_t dimensionIn();

    virtual size_t dimensionOut();

    virtual bool isSelfAdjoint();

    virtual Matrix callAdjoint(Matrix& x);






private:



    LinearOperator& A;
    LinearOperator& B;
};

#endif	/* OPCOMPOSITION_H */

