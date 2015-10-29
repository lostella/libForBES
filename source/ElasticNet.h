/* 
 * File:   ElasticNet.h
 * Author: Pantelis Sopasakis
 *
 * Created on October 28, 2015, 7:43 PM
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

#ifndef ELASTICNET_H
#define	ELASTICNET_H

#include "Function.h"
#include "FunctionOntologicalClass.h"


/**
 * \class ElasticNet
 * \brief Elastic net regularization function
 * \version version 0.1
 * \ingroup Functions
 * \date Created on October 28, 2015, 7:43 PM
 * \author Pantelis Sopasakis
 * 
 * This class extends the class Function and implements the elastic net regularization
 * function which is a function \f$g:\mathbb{R}^n\to\mathbb{R}_{+}\f$ defined as 
 * 
 * \f[
 *  g(x) = \mu \|x\|_1 + \frac{\lambda}{2}\|x\|_2^2.
 * \f]
 * 
 * Its proximal operator is then defined as
 * 
 * \f[
 * \mathrm{prox}_{\gamma g}(v) = \frac{t(v)}{1+\lambda \gamma}.
 * \f]
 * 
 * In order to define function \f$t::\mathbb{R}^n\to\mathbb{R}\f$ we need first to introduce 
 * the following function \f$\psi:\mathbb{R}^n\times\mathbb{R}^n \to \mathbb{R}^n\f$ as
 * 
 * \f[
 *  (\psi(z,v))_i = \begin{cases}
 * z_1,&\text{ if } v_i>0\\
 * -z_i,&\text{ if } v_i>0\\
 * 0,&\text{ otherwise}
 * \end{cases}
 * \f]
 * 
 * We also denote the element-wise absolute value of a vector \f$x\in\mathbb{R}^n\f$ by \f$|x|\f$ and
 * we define the function 
 * 
 * \f[
 *  \tau(v) = \max(0, |v|-\gamma\mu).
 * \f]
 * 
 * Then, function \f$t\f$ is given by
 * 
 * \f[
 * t(v) = \psi(\tau(v), v).
 * \f]
 * 
 * This class also allows the computation of the value of the function at a proximal
 * point, that is \f$g(\mathrm{prox}_{\gamma g}(v))\f$, as a function of \f$v\f$.
 * 
 * This is given by
 * 
 * \f[
 *  g(\mathrm{prox}_{\gamma g}(v)) = \mu 1' \frac{t(v)}{1+\lambda \gamma} + \frac{\lambda}{2}\|\frac{t(v)}{1+\lambda \gamma}\|_2^2.
 * \f]
 * 
 * This a simple MATLAB implementation of this function:
 * 
 * \code
 * function [prox, g] = elasticNetProx(x, gam, mu, lam)
 * uz = max(0, abs(x)-gam*mu)/(1+lam*gam);
 * prox = sign(x).*uz;
 * g = mu*sum(uz)+(0.5*lam)*(uz'*uz);
 * \endcode
 */
class ElasticNet : public Function{
public:
    
    ElasticNet(double lambda, double mu);
    
    virtual ~ElasticNet();

    virtual int call(Matrix& x, double& f);

    virtual int callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox);
    
    virtual int callProx(const Matrix& x, double gamma, Matrix& prox);
    



private:
    double mu;
    double lambda;

};

#endif	/* ELASTICNET_H */

