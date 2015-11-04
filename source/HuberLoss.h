/* 
 * File:   HuberLoss.h
 * Author: Pantelis Sopasakis
 *
 * Created on October 30, 2015, 1:57 AM
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

#ifndef HUBERLOSS_H
#define	HUBERLOSS_H

#include "Function.h"
#include <math.h>

/**
 * \class HuberLoss
 * \brief Huber loss function
 * \version version 0.1
 * \ingroup Functions
 * \date Created on October 30, 2015, 1:57 AM
 * \author Pantelis Sopasakis
 * 
 * This class extends the class Function and implements the Huber loss
 * function which is a function \f$f:\mathbb{R}^n\to\mathbb{R}\f$ defined as 
 * 
 * \f[
 *  f(x) = \sum_{i=1}^{n} l_{\delta}(x_i),
 * \f]
 *
 * where
 * 
 * \f[
 *  l_{\delta}(x_i) = \begin{cases}
 *   \frac{1}{2\delta}\|x_i\|^2,&\text{if } |x_i|\leq \delta \\
 *   |x_i| - \frac{\delta}{2},&\text{otherwise}
 * \end{cases}
 * \f]
 * 
 * The gradient of the Huber loss function is given by
 * 
 * \f[
 *  \nabla f(x)_i = \begin{cases}
 *   \frac{|x_i|}{\delta},&\text{if } |x_i|\leq \delta \\
 *   \mathrm{sign}(x_i),&\text{otherwise}
 * \end{cases}
 * \f]
 * 
 * 
 */
class HuberLoss : public Function {
public:
    
    /**
     * Create a new instance of the Huber loss function with given parameter delta
     * @param delta parameter delta
     */
    explicit HuberLoss(double delta);
    
    /**
     * Destructor
     */
    virtual ~HuberLoss();
    
    virtual int call(Matrix& x, double& f);
    
    virtual int call(Matrix& x, double& f, Matrix& grad);
    
    virtual FunctionOntologicalClass category();

private:
    double m_delta;

};

#endif	/* HUBERLOSS_H */

