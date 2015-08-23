/* 
 * File:   IndBox.h
 * Author: chung
 *
 * Created on July 26, 2015, 5:22 PM
 */

#ifndef INDBOX_H
#define	INDBOX_H

#include "Matrix.h"
#include "Function.h"

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

