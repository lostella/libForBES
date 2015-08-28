/*
 * File:   IndBox.h
 * Author: Pantelis Sopasakis
 *
 * Created on July 26, 2015, 5:22 PM
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

#ifndef INDBOX_H
#define	INDBOX_H

#include "Matrix.h"
#include "Function.h"

/**
 * \class IndBox
 * \brief Indicator of a box
 * \version version 0.1
 * \ingroup Functions
 * \date Created on July 26, 2015, 5:22 PM
 * \author Pantelis Sopasakis
 * 
 * This class extends the class Function and implements the indicator function of
 * a box. 
 * 
 * A box is a set of the form \f$B_{[l,u]}=\{x \in \mathbb{R}^n: l\leq x \leq u\}\f$
 * and its indicator is defined as 
 * \f[
 * \delta (x\mid B_{[l,u]}) = \begin{cases}
 * 0, \text{ if } x\in B_{[l,u]},\\
 * \infty, \text{ otherwise}
 * \end{cases}
 * \f]
 * 
 * The proximal operator of this function is defined as 
 * \f[
 * \delta^*(x^*|B_{[l,u]}) = \mathrm{mid}(x^*; l, u) = \min (\max(x^*, l), u)
 * \f]
 * These two functions are implemented and are available through #call and 
 * #callProx.
 */
class IndBox : public Function {
public:

    IndBox(double& uniform_lb, double& uniform_ub);

    IndBox(Matrix& lb, Matrix& ub);

    virtual ~IndBox();
    
    virtual int call(Matrix& x, double& f);
    
    virtual int category();
    
    virtual int callProx(const Matrix& x, double gamma, Matrix& prox, double& f_at_prox);

    virtual int callProx(const Matrix& x, double gamma, Matrix& prox);


protected:

    void SetLb(Matrix* lb);
    
    void SetUb(Matrix* ub);

    void SetUniform_lb(double* uniform_lb);

    void SetUniform_ub(double* uniform_ub);

    Matrix* GetLb() const;

    Matrix* GetUb() const;

    double* GetUniform_lb() const;

    double* GetUniform_ub() const;

private:

    Matrix* m_lb = NULL;
    Matrix* m_ub = NULL;
    double* m_uniform_lb = NULL;
    double* m_uniform_ub = NULL;


};

#endif	/* INDBOX_H */

