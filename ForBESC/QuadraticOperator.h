/* 
 * File:   QuadraticOperator.h
 * Author: Pantelis Sopasakis
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
 * A quadratic operator is a function of the form \f$Q(x) = \frac{1}{2} x' T(x)\f$,
 * where \f$T\f$ is a given linear operator. 
 * 
 * Operator \f$T:\mathbb{R}^n \to \mathbb{R}^n\f$ produces 
 * a column vector of the same dimension as \f$x\f$. 
 * 
 * It is assumed that 
 * \f$T\f$ is a symmetric operator, therefore the derivative of \f$Q(x)\f$
 * is computed as \f$Q'(x) = T(x)\f$.
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

