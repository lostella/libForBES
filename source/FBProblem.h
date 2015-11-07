#ifndef FBPROBLEM_H
#define	FBPROBLEM_H

#include "Function.h"
#include "LinearOperator.h"

/**
 * \class FBProblem
 * \brief Low level FB problem description
 * \version version 0.0
 * \ingroup FBSolver-group
 * 
 * FBProblem is a low-level call which provides access to the specifications of a 
 * generic optimization problem in the following form:
 * 
 * \f[
 *  \min_{x\in X} F(x),
 * \f]
 * 
 * where \f$X\f$ is a linear space, either \f$\mathbb{R}^n\f$ or \f$\mathbb{R}^{n\times m}\f$
 * 
 */
class FBProblem {
private:
    Function * m_f1;
    Function * m_f2;
    Function * m_g;
    LinearOperator * m_L1;
    LinearOperator * m_L2;
    Matrix * m_d1;
    Matrix * m_d2;
    Matrix * m_lin;

public:

    FBProblem(
            Function& fun_f1,
            LinearOperator& L_1,
            Matrix& d_1,
            Function& fun_f2,
            LinearOperator& L_2,
            Matrix& d_2,
            Matrix& linear,
            Function& fun_g);

    Function * f1();

    LinearOperator * L1();

    Matrix * d1();

    Function * f2();

    LinearOperator * L2();

    Matrix * d2();

    Matrix * lin();

    Function * g();


    virtual ~FBProblem();

};


#endif /* FBPROBLEM_H */
