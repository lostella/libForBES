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

class QuadraticOperator : public Function {
public:

    QuadraticOperator(LinearOperator& T) :
    Function(), T(T) {
    }

    QuadraticOperator(const QuadraticOperator& other) :
    Function(other), T(other.T) {
    }

    virtual ~QuadraticOperator();
    

    virtual int call(Matrix& x, double& f);

    virtual int category();


private:
    LinearOperator& T;

};

#endif	/* QUADRATICOPERATOR_H */

