/* 
 * File:   HingeLoss.h
 * Author: chung
 *
 * Created on October 29, 2015, 10:49 PM
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

#ifndef HINGELOSS_H
#define	HINGELOSS_H

#include "Function.h"

/**
 * \class HingeLoss
 * \brief %Hinge loss function
 * \version version 0.1
 * \ingroup Functions
 * \date Created on October 29, 2015, 10:49 PM
 * \author Pantelis Sopasakis
 * 
 * This class extends the class Function and implements the Hinge loss
 * function which is a function \f$f:\mathbb{R}^n\to\mathbb{R}\f$ defined as 
 * 
 * \f[
 *  f(x) = \mu \sum_{i=1}^{n} \max(0, 1-Bx),
 * \f]
 *
 * where \f$B=\mathrm{diag}(b)\f$ is a diagonal matrix - where \f$b\in\mathbb{R}^n\f$.
 * 
 * Note that in the above expression, \f$z=Bx\in\mathbb{R}^n\f$ is a vector and \f$1-z\f$ 
 * is the vector whose i-th entry is \f$1-z_i\f$. What is more, \f$\max(0, 1-Bx)\f$
 * is also a vector whose i-th entry is \f$\max(0, 1-b_ix_i)\f$.
 * 
 * The proximal operator of the Hinge loss function is given by
 * 
 * \f[
 *  \mathrm{prox}_{\gamma f}(v)_i = 
 *   \begin{cases}
 *     b_i \min(1, b_i + \gamma \mu), &\text{if } b_i x_i < 1\\
 *     x_i, &\text{otherwise}
 *   \end{cases}
 * \f]
 * 
 */
class HingeLoss : public Function {
    
public:

    /**
     * Create a new instance of HingleLoss providing the parameters b and mu.
     * @param b Vector b
     * @param mu Parameter \f$ \mu\f$
     */
    HingeLoss(Matrix& b, double mu);

    /**
     * Create a new instance of HingleLoss providing the parameter b while it is
     * assumed that mu is equal to 1.
     * @param b Vector b
     */
    HingeLoss(Matrix& b);

    /**
     * Destructor
     */
    virtual ~HingeLoss();

    virtual int call(Matrix& x, double& f);
    
    virtual int callProx(const Matrix& x, double gamma, Matrix& prox);
    
    virtual int callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    virtual FunctionOntologicalClass category();


private:

    double m_mu;
    Matrix * m_b;

};

#endif	/* HINGELOSS_H */

