/* 
 * File:   SumOfNorm2.h
 * Author: Pantelis Sopasakis
 *
 * Created on November 6, 2015, 1:23 AM
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

#ifndef SUMOFNORM2_H
#define	SUMOFNORM2_H

#include "Function.h"
#include "Norm.h"



/**
 * \class SumOfNorm2
 * \brief Sum of 2-norms
 * \version version 0.1
 * \ingroup Functions
 * \date Created on November 6, 2015, 1:23 AM
 * \author Pantelis Sopasakis
 * 
 * Let \f$X \cong \mathbb{R}^n\f$ be a vector space and let \f$k\in\mathbb{N}\f$
 * be an integer which divides \f$n\f$ exactly. 
 * 
 * Let us denote \f$m=n/k\f$. Then,
 * 
 * \f[
 * X = \prod_{i=1}^{m}X_i,
 * \f]
 * 
 * where \f$X_i\cong \mathbb{R}^{k}\f$, all equipped with the norm \f$\mu \|\cdot\|_2\f$
 * (for some constant scalar \f$\mu>0\f$).
 * 
 * Then, space \f$X\f$ is equipped with the norm
 * \f[
 * \|x\| = \sum_{i=1}^{k} \mu\|x_{(i)}\|_2.
 * \f]
 * 
 * This is indeed a norm in \f$X\f$ and makes the space complete.
 * 
 * There is also an inner product we may associate with space \f$X\f$, that is
 * 
 * \f[
 *  \langle x,y\rangle = \mu \sum_{i=1}^{k}  x_{(i)}^{\top}y_{(i)}.
 * \f]
 * 
 * The above norm, known as <em>sum of Euclidean norms</em> defines the dual norm
 * 
 * \f[
 * \|x\|_* = \frac{1}{\mu} \max \{\|x_{(1)}\|_2, \|x_{(2)}\|_2, \ldots, \|x_{(m)}\|_2\},
 * \f]
 * 
 * where \f$x_{(i)}\in\mathbb{R}^k\f$ is the part of \f$x\f$ including its components
 * from \f$(i-1)k+1\f$ to \f$ik-1\f$.
 * 
 * \note
 * This implementation assumes an evenly distributed partitioning of \f$X\f$. This
 * may be restrictive in many applications and the class SeparableSum should be
 * used instead.
 * 
 * \note
 * This function is used in group LASSO problems.
 * 
 */ 
class SumOfNorm2 : public Norm {
public:
    using Function::call;
    
    /**
     * Defines a sum-of-norms function with a given partition length \f$k\f$
     * @param k partition length
     */
    explicit SumOfNorm2(size_t k);
    
    /**
     * Defines a sum-of-norms function with a given partition length \f$k\f$ and
     * a given scaling parameter
     * 
     * @param mu scaling parameter
     * @param k partition length
     */
    SumOfNorm2(double mu, size_t k);

    /**
     * Default destructor. This will <code>delete</code> an internal private member
     * <code>m_norm</code> which is an instance of <code>Function</code> used to 
     * compute the norm-2 values of the subvectors.
     */
    virtual ~SumOfNorm2();
    
    virtual int call(Matrix& x, double& f);

    /**
     * The proximal operator \f$\mathrm{prox}_{\gamma f}\f$ for the sum-of-norms
     * function is computed as follows
     * 
     * \f[
     *  (\mathrm{prox}_{\gamma f}(v))_{(i)} = \left(1-\frac{\gamma}{\|v_{(i)}\|_2}\right)v_{(i)}
     * \f]
     * 
     * @param v Vector \f$v\f$ where the proximal operator is evaluated
     * @param gamma Parameter \f$\gamma\f$ of the proximal operator 
     * @param prox The result of this operation
     * 
     * @return 
     * status code which is equal to <code>STATUS_OK=0</code> if the computation
     * has succeeded without any problems, <code>STATUS_UNDEFINED_FUNCTION=2</code> if
     * this function is not defined by the derived class and <code>STATUS_NUMERICAL_PROBLEMS=1</code>
     * if some numerical problems prevented the computation of a reliable result. 
     * Custom implementations are allowed to return other non-zero error/warning
     * status codes.
     * 
     * \exception std::invalid_argument an <code>invalid_argument</code> exception
     * is thrown in case the function argument <code>x</code> and/or <code>prox</code>
     * are of incompatible dimensions.
     * 
     */
    virtual int callProx(Matrix& v, double gamma, Matrix& prox);

    virtual int callProx(Matrix& v, double gamma, Matrix& prox, double& f_at_prox);    

    /**
     * Computes the dual norm of sum-of-norms
     * @param x
     * @param norm
     * @return 
     */
    virtual int dualNorm(Matrix& x, double& norm);

    virtual FunctionOntologicalClass category();

private:
    /**
     * Parameter mu (scaling)
     */
    double m_mu;
    
    /**
     * The size of each chung of <code>x</code>.
     */
    size_t m_partition_length;
    
    /**
     * Pointer to a Norm2 function.
     */
    Function * m_norm2;

};

#endif	/* SUMOFNORM2_H */

