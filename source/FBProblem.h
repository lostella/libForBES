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
 *  \min_{x\in X} f_1(L_1 x + d_1) + f_2(L_2 x + d_2) + \langle l,x \rangle + g(x),
 * \f]
 * 
 * where \f$X\f$ is either \f$\mathbb{R}^n\f$ or \f$\mathbb{R}^{n\times m}\f$
 * and \f$f_1\f$ is a smooth quadratic function, \f$f_2\f$ is a smooth function,
 * \f$L_1, L_2\f$ are linear operators and \f$d_1, d_2\f$ are elements of the
 * respective target spaces, \f$l\in X\f$ and \f$g\f$ is a closed, proper, convex
 * function over X.
 * 
 * 
 * \example qp_box.cpp
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

    void init();
    

public:

    /**
     * Allocates an FBProblem given all the details. 
     * 
     * Defines the optimization problem
     * 
     * \f[
     *  \mathrm{minimize}\  f_1(L_1 x + d_1) + f_2(L_2 x + d_2) + \langle l,x \rangle + g(x),
     * \f]
     * 
     * @param fun_f1 %Quadratic function \f$f_1\f$
     * @param L_1 Linear operator \f$L_1\f$
     * @param d_1 Constant matrix/vector \f$d_1\f$
     * @param fun_f2 Smooth %function \f$f_2\f$
     * @param L_2 Linear operator \f$L_2\f$
     * @param d_2 Constant matrix/vector \f$d_2\f$ 
     * @param linear Constant \f$l\f$
     * @param fun_g Closed proper convex function \f$g\f$
     *
     * \todo Test for "quadraticness" of fun_f1 and fun_f2, and appropriately assign
     * internally all the input arguments internally. If, for example, they are both
     * quadratic or both non-quadratic, then they should be summed.
     */
    FBProblem(
            Function& fun_f1,
            LinearOperator& L_1,
            Matrix& d_1,
            Function& fun_f2,
            LinearOperator& L_2,
            Matrix& d_2,
            Matrix& linear,
            Function& fun_g);
    
    /**
     * Allocates an FBProblem given only one smooth term and the affine
     * map composed with it. 
     * 
     * Defines the optimization problem
     * 
     * \f[
     *  \mathrm{minimize}\  f(Lx + d) + g(x)
     * \f]
     *
     * \todo Test for "quadraticness" of fun_f, and appropriately assign internally
     * all the input arguments.
     * 
     * @param fun_f Function \f$f\f$
     * @param L linear operator \f$L\f$
     * @param d constant matrix/vector \f$d\f$
     * @param fun_g Closed proper convex function \f$g\f$
     */
    FBProblem(
            Function& fun_f,
            LinearOperator& L,
            Matrix& d,
            Function& fun_g);

    /**
     * Allocates an FBProblem given only one smooth term and the linear map
     * composed with it.
     * 
     * Defines the optimization problem
     * 
     * \f[
     *  \mathrm{minimize}\  f(Lx) + g(x)
     * \f]
     * 
     * @param fun_f Function \f$f\f$
     * @param L Linear operator \f$L\f$
     * @param fun_g Closed proper convex function \f$g\f$
     *
     * \todo Test for "quadraticness" of fun_f, and appropriately assign internally
     * all the input arguments.
     */
    FBProblem(
            Function& fun_f,
            LinearOperator& L,
            Function& fun_g);

    /**
     * Allocates an FBProblem given only one smooth term.
     * 
     * \f[
     *  \mathrm{minimize}\ f(x) + g(x)
     * \f]
     * 
     * @param fun_f Function \f$f\f$
     * @param fun_g Closed proper convex function \f$g\f$
     *
     * \todo Test for "quadraticness" of fun_f, and appropriately assign internally
     * all the input arguments.
     */
    FBProblem(
            Function& fun_f,
            Function& fun_g);

    /**
     * Quadratic function in the problem objective.
     *
     * @return Pointer to the quadratic Function in the problem.
     */
    Function * f1();

    /**
     * Non-quadratic function in the problem objective.
     *
     * @return Pointer to the non-quadratic Function in the problem.
     */
    Function * f2();

    /**
     * Linear operator composed with the quadratic function.
     *
     * @return Pointer to the linear operator.
     */
    LinearOperator * L1();

    /**
     * Linear operator composed with the non-quadratic function.
     *
     * @return Pointer to the linear operator.
     */
    LinearOperator * L2();

    /**
     * Affine term composed with the quadratic function.
     *
     * @return Pointer to Matrix.
     */
    Matrix * d1();

    /**
     * Affine term composed with the non-quadratic function.
     *
     * @return Pointer to Matrix.
     */
    Matrix * d2();

    /**
     * Linear term in the cost.
     *
     * @return Pointer to Matrix.
     */
    Matrix * lin();

    /**
     * The proper, closed, convex function g in the cost.
     *
     * @return Pointer to Function.
     */
    Function * g();

    virtual ~FBProblem();

};


#endif /* FBPROBLEM_H */
