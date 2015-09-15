/* 
 * File:   OpAdjoint.h
 * Author: Pantelis Sopasakis
 *
 * Created on September 15, 2015, 2:57 PM
 */

#ifndef OPADJOINT_H
#define	OPADJOINT_H

#include "LinearOperator.h"

/**
 * \class OpAdjoint
 * \brief Adjoint of a given linear operator
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on September 15, 2015, 2:57 PM
 * 
 * \ingroup LinOp
 */
class OpAdjoint : public LinearOperator {
public:

    OpAdjoint(LinearOperator& op);

    virtual ~OpAdjoint();

    virtual Matrix call(Matrix& x);

    virtual Matrix callAdjoint(Matrix& x);

    virtual size_t dimensionIn();

    virtual size_t dimensionOut();

    virtual bool isSelfAdjoint();

private:

    LinearOperator& m_originalOperator;

};

#endif	/* OPADJOINT_H */

