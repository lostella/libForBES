/* 
 * File:   SeparableSum.h
 * Author: Pantelis Sopasakis
 *
 * Created on October 30, 2015, 7:22 PM
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

#ifndef SEPARABLESUM_H
#define	SEPARABLESUM_H

#include <vector>
#include <map>

#include "Function.h"

/**
 * \class SeparableSum
 * \brief Separable sum
 * \version version 0.1
 * \ingroup Functions
 * \date Created on October 30, 2015, 7:22 PM
 * \author Pantelis Sopasakis
 * 
 * For \f$x\in\mathbb{R}^n\f$ 
 * and \f$\mathcal{I}\f$ be a set of unique indices in \f$\mathbb{N}_{[1,n]}\f$ given by
 * \f$\mathcal{I}=\{i^1, i^2, \ldots, i^r\}\f$, define \f$x_{\mathcal{I}} = (x_{i^1},\ldots, x_{i^r})\f$.
 * 
 * Take a function \f$f:\mathbb{R}^n\to\mathbb{R}\f$ and sets of indices 
 * \f$\mathcal{I}_1, \ldots, \mathcal{I}_2\f$ so that \f$\mathcal{I}_i \cap \mathcal{I}_j = \emptyset\f$
 * whenever \f$i \neq j\f$ so that
 * 
 * \f[
 *  f(x) = \sum_{i=1}^{s} f_i(x_{\mathcal{I}_j}).
 * \f]
 * 
 * Such a function is called a separable sum.
 * 
 * The proximal operator of a separable sum can be computed as follows:
 * 
 * \f[
 * (\mathrm{prox}_{\gamma f}(v))_{\mathcal{I}_j} = \mathrm{prox}_{\gamma f_j}(v_{\mathcal{I}_j}).
 * \f]
 * 
 * As a result, it suffices to be able to compute \f$ \mathrm{prox}_{\gamma f_j}(\cdot)\f$
 * to compute the proximal of a separable sum.
 */
class SeparableSum : public Function {
public:    
    using Function::call;
    using Function::callConj;

    virtual ~SeparableSum();

    /**
     * Constructor for an instance of a separable sum function given its characteristic
     * map which maps Function pointers to sets of indices. In particular, to each 
     * function \f$f_{i}\f$ there is an associated set of indices \f$\mathcal{I}_i\f$
     * so that \f$f_{i}\f$ is a function of \f$x_{\mathcal{I}_i}\f$.
     * 
     * @param fun_idx_map map of Functions to sets of indices.
     */
    explicit SeparableSum(std::map<Function*, std::vector<size_t>*> fun_idx_map);

    virtual int call(Matrix& x, double& f);

    virtual int call(Matrix& x, double& f, Matrix& grad);
    
    /**
     * Computes the conjugate of a separable sum function using the property that for
     * a separable sum of the form
     * 
     * \f[
     *  f(x) = \sum_{i=1}^{s} f_i(x_{\mathcal{I}_j}).
     * \f]
     * 
     * the conjugate is given by
     * 
     * \f[
     *  f^*(y) = \sum_{i=1}^{s} f_i^*(y_{\mathcal{I}_j}).
     * \f]
     * 
     * @param x given vector of appropriate dimensions
     * @param f_star value of the conjugate function to be computed
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     */
    virtual int callConj(Matrix& x, double& f_star);


    /**
     * Computes the proximal operator of a separable sum of functions using the
     * formula
     * 
     * \f[
     * (\mathrm{prox}_{\gamma f}(v))_{\mathcal{I}_j} = \mathrm{prox}_{\gamma f_j}(v_{\mathcal{I}_j}).
     * \f]
     * 
     * @param x given vector of appropriate dimensions
     * @param gamma gamma parameter of the proximal operator
     * @param prox the result of the prox operator
     * @return 
     * status code which is equal to <code>STATUS_OK</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     * 
     * \exception std::invalid_argument if <code>x</code> is not a column vector
     * or if it has incompatible dimensions.
     * 
     */
    virtual int callProx(Matrix& x, double gamma, Matrix& prox);

    /**
     * Computes the proximal operator of a separable sum of functions using the
     * formula
     * 
     * \f[
     * (\mathrm{prox}_{\gamma f}(v))_{\mathcal{I}_j} = \mathrm{prox}_{\gamma f_j}(v_{\mathcal{I}_j}).
     * \f]
     * 
     * This method computes both \f$\mathrm{prox}_{\gamma f}(v)\f$ and \f$f(\mathrm{prox}_{\gamma f}(v))\f$.
     * 
     * @param x given vector of appropriate dimensions
     * @param gamma gamma parameter of the proximal operator
     * @param prox the result of the prox operator
     * @param f_at_prox the value of f at the proximal point
     * @return 
     * status code which is equal to <code>STATUS_OK</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     * 
     * \exception std::invalid_argument if <code>x</code> is not a column vector
     * or if it has incompatible dimensions.
     * 
     */
    virtual int callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox);
   
    virtual FunctionOntologicalClass category();


private:

    std::map<Function*, std::vector<size_t> * > m_fun_idx_map;

};

#endif	/* SEPARABLESUM_H */

