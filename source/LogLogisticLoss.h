/* 
 * File:   LogLogisticLoss.h
 * Author: Pantelis Sopasakis
 *
 * Created on October 29, 2015, 5:08 PM
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

#ifndef LOGLOGISTICLOSS_H
#define	LOGLOGISTICLOSS_H

#include "Function.h"
#include <math.h>

/**
 * \class LogLogisticLoss
 * \brief Log-logistic loss function
 * \version version 0.1
 * \ingroup Functions
 * \date Created on October 29, 2015, 5:08 PM
 * \author Pantelis Sopasakis
 * 
 * This class extends the class Function and implements the log-logistic loss
 * function which is a function \f$f:\mathbb{R}^n\to\mathbb{R}\f$ defined as 
 * 
 * \f[
 *  f(x) = \mu \sum_{i=1}^{n} -\ln \sigma(x_i),
 * \f]
 * 
 * where
 * 
 * \f[
 *  \sigma(z) = \frac{e^{z}}{1+e^{z}}.
 * \f]
 * 
 * The gradient of \f$f\f$, \f$\nabla f:\mathbb{R}^n \to \mathbb{R}^n \f$ is computed element-wise by
 * 
 * \f[
 * \nabla f(x)_i = \mu ( \sigma(x_i) - 1 ).
 * \f]
 * 
 * The Hessian of this function is a diagonal matrix with diagonal elements
 * 
 * \f[
 *  (\nabla^2 f(x))_{ii} = \frac{\mu e^{x_i}}{(1+e^{x_i})^2}.
 * \f]
 */

class LogLogisticLoss : public Function {
public:
    
    /**
     * Create an instance of LogLogisticLoss assuming \f$\mu=1\f$
     */
    LogLogisticLoss(); // with mu = 1    
    
    /**
     * Create an instance of LogLogisticLoss using a given value for \f$\mu\f$
     * @param mu Parameter \f$\mu\f$ (positive)
     */
    explicit LogLogisticLoss(double mu);

    virtual ~LogLogisticLoss();
    
    virtual int call(Matrix& x, double& f, Matrix& grad);

    virtual int call(Matrix& x, double& f);

    virtual FunctionOntologicalClass category();



private:
    double mu;

};

#endif	/* LOGLOGISTICLOSS_H */

