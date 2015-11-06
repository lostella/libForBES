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

    int m_n;
    int m_m1;
    int m_m2;

public:

    FBProblem(
            int n,
            Function * f1,
            LinearOperator * L1,
            Matrix * d1,
            Function * f2,
            LinearOperator * L2,
            Matrix * d2,
            Matrix * lin,
            Function * g);

    Function * getf1() {
        return m_f1;
    }

    LinearOperator * getL1() {
        return m_L1;
    }

    Matrix * getd1() {
        return m_d1;
    }

    Function * getf2() {
        return m_f2;
    }

    LinearOperator * getL2() {
        return m_L2;
    }

    Matrix * getd2() {
        return m_d2;
    }

    Matrix * getlin() {
        return m_lin;
    }

    Function * getg() {
        return m_g;
    }

    int getn() {
        return m_n;
    }

    int getm1() {
        return m_m1;
    }

    int getm2() {
        return m_m2;
    }

    virtual ~FBProblem() {
    }

};


#endif /* FBPROBLEM_H */
