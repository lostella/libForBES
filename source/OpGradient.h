/* 
 * File:   OpGradient.h
 * Author: chung
 *
 * Created on September 16, 2015, 1:35 AM
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

#ifndef OPGRADIENT_H
#define	OPGRADIENT_H

#include "LinearOperator.h"

/**
 * \class OpGradient
 * \brief Discrete gradient transform
 * \version 0.0-tentative
 * \author Pantelis Sopasakis
 * \date Created on September 16, 2015, 1:35 AM
 * 
 * \ingroup LinOp
 * 
 * For a vector \f$x=(x_0, x_1, \ldots, x_{n-1})\in\mathbb{R}^n\f$ with \f$n\geq 2\f$, 
 * we define its <em>gradient</em> as the vector \f$G(x)\in\mathbb{R}^{n-1}\f$ with 
 * \f[
 * (G(x))_{k} = x_{k+1}-x_{k}, 
 * \f]
 * for \f$k=0,\ldots, n-2\f$. Operator \f$G:\mathbb{R}^{n}\to\mathbb{R}^{n-1}\f$ is a linear transform. 
 * The <em>total variation</em> of \f$x\f$ is then defined as
 * \f[
 * \mathrm{TV}(x) = \|G(x)\|_1.
 * \f]
 * The conjugate operator \f$G^*:\mathbb{R}^{n-1}\to \mathbb{R}^n\f$, which is referred 
 * to as the <em>divergence</em> of \f$x\f$, is such that
 * \f[
 * (G^*(y))_0=-y_0,
 * \f]
 * and for \f$k=1,\ldots, n-2\f$ it is
 * \f[
 * (G^*(y))_k=y_k-y_{k+1},
 * \f]
 * and 
 * \f[
 * (G^*(y))_{n-1}=y_{n-2}.
 * \f]
 * One can verify that for \f$x\in\mathbb{R}^n\f$ and \f$y\in\mathbb{R}^{n-1}\f$ we have
 * \f$\langle y, G(x) \rangle = \langle G^*(y), x \rangle\f$.
 */
class OpGradient : public LinearOperator {
public:
    OpGradient();

    OpGradient(size_t n);

    virtual ~OpGradient();

    virtual Matrix call(Matrix& x);

    virtual Matrix callAdjoint(Matrix& x);

    virtual size_t dimensionIn();

    virtual size_t dimensionOut();

    virtual bool isSelfAdjoint();


private:

    size_t m_dimension;

};

#endif	/* OPGRADIENT_H */

