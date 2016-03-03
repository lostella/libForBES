/* 
 * File:   LQCost.h
 * Author: chung
 *
 * Created on March 3, 2016, 2:00 PM
 */

#ifndef LQCOST_H
#define	LQCOST_H

#include "Function.h"
#include "Matrix.h"
#include "CholeskyFactorization.h"

/**
 * \class LQCost
 * \brief Smooth part of a linear-quadratic optimal control problem
 * \version version 0.1
 * \ingroup Functions
 * \date Created on March 3, 2016, 2:00 PM
 * \author Pantelis Sopasakis
 * 
 */
class LQCost : public Function {
public:

    LQCost();;

    virtual ~LQCost();

    virtual int callConj(Matrix& x, double& f_star, Matrix& grad);

    virtual int callConj(Matrix& x, double& f_star);

    virtual FunctionOntologicalClass category();


private:

    /* Problem data */
    Matrix * m_A;
    Matrix * m_B;
    Matrix * m_f;
    Matrix * m_Q;
    Matrix * m_QN;
    Matrix * m_R;
    Matrix * m_S;
    Matrix * m_r;
    Matrix * m_q;
    Matrix * m_qN;
    size_t m_N;

    /* Current state */
    Matrix * m_p;

    /* Internal data (output of factor step) */
    Matrix * m_L;
    Matrix * m_K;
    Matrix * m_Rbar;
    Matrix * m_d;
    Matrix * m_s;
    CholeskyFactorization * m_RbarFactor;


    /**
     * Updates the data of the factor step. This is to be invoked once. 
     * @return LibForBES status code 
     */
    int factor_step();

    /**
     * Initializes everything to \c NULL.
     */
    void nullify_all();



};

#endif	/* LQCOST_H */

