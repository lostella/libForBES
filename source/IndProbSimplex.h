/* 
 * File:   IndProbSimplex.h
 * Author: chung
 *
 * Created on January 15, 2016, 2:37 AM
 */

#ifndef INDPROBSIMPLEX_H
#define	INDPROBSIMPLEX_H

#include "Function.h"

/**
 * \class IndProbSimplex
 * \brief %Indicator of the probability simplex
 * \version 0.1
 * \author Pantelis Sopasakis
 * \date Created on January 15, 2016, 2:37 AM
 * 
 * This class implements the indicator of the probability simplex, that is
 * 
 * \f[
 * f(x) = \delta(x\mid P) = 
 * \begin{cases}
 *  0, &\text{if } x\in P\\
 * +\infty,&\text{otherwise}
 * \end{cases}
 * \f]
 * 
 * where \f$P\f$ is the set
 * 
 * \f[
 * P = \{x\in\mathbb{R}^n: x\geq 0, 1'x = 1\}.
 * \f]
 * 
 * The projection on this set is given by
 * 
 * \f[
 * \mathrm{proj}(x\mid P) = (x-t1)_+,
 * \f]
 * 
 * where \f$(z)_+=\max \{z, 0\}\f$ (element-wise) and \f$\max_i x_i - 1 \leq t \leq \max_i x_i\f$ 
 * is a scalar such that
 * 
 * \f[
 * 1'(x-t1)_+ = 1 \Leftrightarrow \sum_{i=1}^{n}(x_i-t)_+ = 1.
 * \f]
 * 
 * This can be determined by bisection.
 */
class IndProbSimplex : public Function {
public:
    using Function::call;    
    
    /**
     * Constructs a new instance of the indicator function of the probability
     * simples.
     */
    IndProbSimplex();

    /**
     * Default destructor.
     */
    virtual ~IndProbSimplex();
    
    virtual int call(Matrix& x, double& f);

    virtual int callProx(Matrix& x, double gamma, Matrix& prox);
    
    virtual int callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    virtual FunctionOntologicalClass category();


private:

};

#endif	/* INDPROBSIMPLEX_H */

