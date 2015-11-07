/* 
 * File:   ConjugateFunction.h
 * Author: chung
 *
 * Created on November 7, 2015, 3:15 AM
 */

#ifndef CONJUGATEFUNCTION_H
#define	CONJUGATEFUNCTION_H

#include "Function.h"

/**
 * \class ConjugateFunction
 * \brief Conjugate of a given function
 * \version version 0.1
 * \ingroup Functions
 * \date Created on November 7, 2015, 3:15 AM
 * \author Pantelis Sopasakis
 * 
 * This class is a utility which provides access to the convex conjugate of a 
 * function which is passed to the constructor of ConjugateFunction by reference.
 * 
 * For a given function \f$s:X\to\mathbb{R}\cup \{+\infty\}\f$, this class 
 * provides access to \f$f(x) = s^*(x)\f$. Then,
 * 
 * 1. %Function \f$f(x) = s^*(x)\f$
 * 2. Gradient \f$\nabla f(x) = \nabla s^*(x)\f$
 * 3. Conjugate function \f$f^*(x) = s(x)\f$
 * 4. Conjugate gradient \f$\nabla f^*(x) = \nabla s(x)\f$
 * 5. Proximal operator \f$\mathrm{prox}_{\gamma s^*}(v) = v - \mathrm{prox}_{\gamma s}(v)\f$
 * 
 */
class ConjugateFunction : public Function {
public:

    explicit ConjugateFunction(const Function& funct);

    virtual ~ConjugateFunction();

    /**
     * For a given function \f$s\f$, provided in the constructor, this method 
     * computes its conjugate \f$s^*\f$.
     * @param x
     * @param f
     * @return 
     */
    virtual int call(Matrix& x, double& f);

    /**
     * For a given function \f$s\f$, provided in the constructor, this method 
     * computes its conjugate \f$s^*\f$ and its gradient \f$\nabla s^*\f$.
     * @param x
     * @param f
     * @return 
     */
    virtual int call(Matrix& x, double& f, Matrix& grad);

    /**
     * For a given function \f$s\f$, provided in the constructor, this method 
     * computes its value \f$s(x)\f$ at a given point \f$x\in X\f$.
     * @param x
     * @param f_star
     * @param grad
     * @return 
     */
    virtual int callConj(const Matrix& x, double& f_star, Matrix& grad);

    virtual int callConj(const Matrix& x, double& f_star);

    /**
     * For a given function \f$s\f$, provided in the constructor, this method 
     * computes the proximal gradient of its conjugate which is given by
     * \f[
     *  \mathrm{prox}_{\gamma s^*}(v) = v - \mathrm{prox}_{\gamma s}(v).
     * \f]
     * @param x
     * @param gamma
     * @param prox
     * @return 
     */
    virtual int callProx(const Matrix& x, double gamma, Matrix& prox);

    /**
     * For a given function \f$s\f$, provided in the constructor, this method 
     * computes the proximal gradient of its conjugate which is given by
     * \f[
     *  \mathrm{prox}_{\gamma s^*}(v) = v - \mathrm{prox}_{\gamma s}(v).
     * \f]
     * and it also computes the value of \f$s^*\f$ at the proximal point, that is
     * \f$s^*(\mathrm{prox}_{\gamma s^*}(v))\f$.
     * 
     * @param x
     * @param gamma
     * @param prox
     * @param f_at_prox
     * @return 
     */
    virtual int callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    virtual FunctionOntologicalClass category();



private:

    Function& m_function;

};

#endif	/* CONJUGATEFUNCTION_H */

