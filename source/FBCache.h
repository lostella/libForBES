/*
 * File:   FBCache.h
 * Author: Lorenzo Stella
 *
 * Created on October 2, 2015
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

#ifndef FBCACHE_H
#define	FBCACHE_H

#include "Matrix.h"
#include "FBProblem.h"

/**
 * \class FBProblem
 * \brief Low level FB problem description
 * \version version 0.0
 * \ingroup FBSolver-group
 * 
 * 
 */
class FBCache {
private:
    int m_flag_evalf;
    int m_flag_gradstep;
    int m_flag_proxgradstep;

    FBProblem & m_prob;
    Matrix & m_x;

    /* Reference to problem details */
    Function * m_f1;
    Function * m_f2;
    Function * m_g;
    LinearOperator * m_L1;
    LinearOperator * m_L2;
    Matrix * m_d1;
    Matrix * m_d2;
    Matrix * m_lin;

    /* Internal storage for computing proximal-gradient steps */
    Matrix * m_y;
    Matrix * m_z;
    Matrix * m_res1x;
    Matrix * m_gradf1x;
    Matrix * m_res2x;
    Matrix * m_gradf2x;
    Matrix * m_gradfx;
    double m_f1x;
    double m_f2x;
    double m_linx;
    double m_fx;
    double m_gz;
    double m_gamma;

    int update_forward_step(double gamma);
    
public:
    FBCache(FBProblem & p, Matrix & myx, double mygamma);

    virtual ~FBCache();

    int update_eval_FBE(double gamma);
    int update_grad_FBE(double gamma);
    int update_eval_f();
    int update_forward_backward_step(double gamma);

};

#endif /* FBCACHE_H */