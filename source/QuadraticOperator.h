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
 * \version 0.1
 * \date Created on July 24, 2015, 8:49 PM
 * 
 *
 * A generalized quadratic function is a function of the form
 * 
 * \f[
 * F(x) = \frac{1}{2}\langle x, T(x) \rangle
 * \f]
 * 
 * where f$T:\mathbb{R}^n \to \mathbb{R}^n\f$ is a linear operator (here, an instance of LinearOperator).
 * This operator is expected to be self-adjoint.
 * 
 * 
 * Here is a very simple example of use where <code>T</code> is an instance of MatrixOperator.
 * 
 * \code
 * double fval;                         // the value of F at x, F(x) - to be computed
 * Matrix grad;                         // the gradient of F at x - to be computed
 * 
 * Matrix Q;
 * Matrix x;
 * 
 * LinearOperator *T = new MatrixOperator(Q);
 * Function *F = new QuadraticOperator(T);
 * 
 * int info = F -> call(x, fval, grad); // compute its value at x, F(x) and its gradient grad(F)(x)
 * 
 * delete T;
 * delete F;
 * \endcode
 *  
 * \ingroup Functions
 * 
 */
class QuadraticOperator : public Function {
public:

    /**
     * Construct a new instance of QuadraticOperator given a linear operator 
     * @param T linear operator
     */
    explicit QuadraticOperator(LinearOperator& T);

    /**
     * Destructor
     */
    virtual ~QuadraticOperator();

    virtual int call(Matrix& x, double& f);
    
    virtual int call(Matrix& x, double& f, Matrix& grad);

    virtual FunctionOntologicalClass category();
    
protected:
    virtual int computeGradient(Matrix& x, Matrix& grad);

private:
    LinearOperator& m_T;

};

#endif	/* QUADRATICOPERATOR_H */

