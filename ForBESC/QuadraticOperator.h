/* 
 * File:   QuadraticOperator.h
 * Author: chung
 *
 * Created on July 24, 2015, 8:49 PM
 */

#ifndef QUADRATICOPERATOR_H
#define	QUADRATICOPERATOR_H

#include "LinearOperator.h"
#include "Function.h"

/**
 * \class QuadraticOperator
 * \brief %Quadratic operator <code>F(x) = 0.5*x'*T(x)</code>, where <code>T</code> is a linear operator
 * \author Pantelis Sopasakis
 * \version 0.0
 * \date Created on July 24, 2015, 8:49 PM
 * 
 * A quadratic operator is a function of the form <code>Q(x) = 0.5*x'*T(x)</code>,
 * where <code>T</code> is a given linear operator. Operator <code>T</code> produces 
 * a column vector of the same dimension as <code>x</code>. It is assumed that 
 * <code>T</code> is a symmetric operator, therefore the derivative of <code>Q(x)</code>
 * is computed as <code>Q'(x) = T(x)</code>.
 * 
 * \ingroup Functions
 * 
 */
class QuadraticOperator : public Function {
public:

    QuadraticOperator(LinearOperator& T) :
    Function(), T(T) {
        if (T.dimensionIn() != T.dimensionOut()) {
            throw std::invalid_argument("T has incompatible dimensions");
        }
    }

    virtual ~QuadraticOperator();

    virtual int call(Matrix& x, double& f);
    virtual int call(Matrix& x, double& f, Matrix& grad);

    virtual int category();

protected:
    virtual int computeGradient(Matrix& x, Matrix& grad);

private:
    LinearOperator& T;

};

#endif	/* QUADRATICOPERATOR_H */

