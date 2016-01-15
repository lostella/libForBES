/* 
 * File:   DistanceToBall2.h
 * Author: chung
 *
 * Created on January 14, 2016, 2:02 AM
 */

#ifndef DISTANCETOBALL2_H
#define	DISTANCETOBALL2_H

#include "Function.h"
#include <cmath>


/**
 * \class DistanceToBall2
 * \brief Distance from the unit ball of R^n
 * \version version 0.1
 * \ingroup Functions
 * \date Created on January 14, 2016, 2:02 AM
 * \author Pantelis Sopasakis
 * 
 * This class extends Function and implements the squared distance to the closed ball
 * of norm-2 in \f$\mathbb{R}^n\f$ of radius \f$\rho\f$, centered at \f$c\in\mathbb{R}^n\f$
 * which is defined by
 * 
 * \f[
 *  B_2(\rho, c) = \{x\in\mathbb{R}^n: \|x-c\|_2 \leq \rho\},
 * \f]
 * 
 * that is, this is the function
 * 
 * \f[
 * f(x) = \frac{w}{2}d(x, B_2(\rho, c))^2 = \frac{w}{2}\|x-\mathrm{proj}(x, B_2(\rho, c))\|^2,
 * \f]
 * 
 * where \f$w>0\f$ is a positive scalar and
 * \f[
 *  \mathrm{proj}(x, B_2) = \mathrm{argmin}\limits_{\|y-c\|\leq \rho}\|y-x\|
 * \f] 
 * is the projection of \f$x\f$ on \f$B_2(\rho, c)\f$. This is computed as
 * 
 * \f[
 * \mathrm{proj}(x, B_2) = \begin{cases}
 * x,   &\text{ if } \|x-c\|\leq \rho\\
 * c+\frac{\rho}{\|x-c\|}(x-c), &\text{ otherwise}
 * \end{cases}
 * \f]
 * 
 * Essentially, the squared distance to ball-2 is given by
 * 
 * \f[
 * f(x) = \begin{cases}
 * 0,  &\text{ if } \|x-c\|\leq \rho,\\
 * \frac{w}{2}(\|x-c\|-\rho)^2,&\text{ otherwise}
 * \end{cases}
 * \f]
 * 
 * 
 * The gradient of \f$f\f$ is given by
 * 
 * \f[
 * \begin{align}
 *  \nabla f (x) &= \frac{w}{2}(x-\mathrm{proj}(x, B_2(\rho, c)))\\
 *               &= \begin{cases}
 *                      0, & \text{if } \|x-c\|\leq \rho,\\
 *                      \left(1-\frac{\rho}{\|x-c\|}\right)(x-c), &\text{otherwise}
 *                  \end{cases}
 * \end{align}
 * \f]
 * 
 */
class DistanceToBall2 : public Function {
public:
    
    using Function::call;
    
    /**
     * Constructs the indicator of the closed norm-2 unit ball centered at the
     * origin, i.e., the function
     * 
     * \f[
     * f(x) = \frac{1}{2}d(x, B_2),
     * \f]
     * 
     * where
     * 
     * \f[
     * B_2 = \{x\in\mathbb{R}^n: \|x\|_2\leq 1\}.
     * \f]
     */
    DistanceToBall2();
    
    
    /**
     * Constructs the indicator of the closed norm-2 unit ball centered at the
     * origin scaled by a given scalar. The function has the form
     * 
     * \f[
     * f(x) = \frac{w}{2}d(x, B_2),
     * \f]
     * 
     * where
     * 
     * \f[
     * B_2 = \{x\in\mathbb{R}^n: \|x\|_2\leq 1\}.
     * \f]
     * 
     * @param w scaling parameter
     */
    explicit DistanceToBall2(double w);
    
    
    /**
     * Constructs the indicator of the closed norm-2 ball centered at the
     * origin with radius \f$\rho\f$ scaled by a given scalar. The function has the form
     * 
     * \f[
     * f(x) = \frac{w}{2}d(x, B_2(\rho)),
     * \f]
     * 
     * where
     * 
     * \f[
     * B_2(\rho) = \{x\in\mathbb{R}^n: \|x\|_2\leq \rho\}.
     * \f]
     * 
     * @param w scaling factor 
     * @param rho radius of the ball
     */
    DistanceToBall2(double w, double rho);
    
    /**
     * Constructs the indicator of the closed norm-2 ball centered at the
     * a point \f$c\in\mathbb{R}^n\f$ with radius \f$\rho\f$ scaled by a 
     * given scalar. The function has the form
     * 
     * \f[
     * f(x) = \frac{w}{2}d(x, B_2(\rho, c)),
     * \f]
     * 
     * where 
     * 
     * \f[
     * B_2(\rho,c) = \{x\in\mathbb{R}^n: \|x-c\|_2\leq \rho\}.
     * \f]
     * 
     * 
     * @param w scaling factor 
     * @param rho radius of the ball
     * @param center center of the ball
     */
    DistanceToBall2(double w, double rho, Matrix &center);
    
    /**
     * Default destructor.
     */
    virtual ~DistanceToBall2();

    virtual int call(Matrix& x, double& f, Matrix& grad);   

    virtual int call(Matrix& x, double& f);

    virtual FunctionOntologicalClass category();
    
private:
    
    /**
     * ball radius
     */
    double m_rho;
    /**
     * scaling factor (the function is scaled by w/2)
     */
    double m_w;
    /**
     * center of the ball (vector).
     * Default is the origin (in which case this pointer is \c NULL).
     */
    Matrix * m_center;

};

#endif	/* DISTANCETOBALL2_H */

