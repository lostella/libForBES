/* 
 * File:   SeparableSum.h
 * Author: Lorenzo Stella
 * Author: Pantelis Sopasakis
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

#include "FBProblem.h"

FBProblem::FBProblem(
        int n,
        Function * f1,
        LinearOperator * L1,
        Matrix * d1,
        Function * f2,
        LinearOperator * L2,
        Matrix * d2,
        Matrix * lin,
        Function * g) {
    this->m_f1 = f1;
    this->m_L1 = L1;
    this->m_d1 = d1;
    this->m_f2 = f2;
    this->m_L2 = L2;
    this->m_d2 = d2;
    this->m_lin = lin;
    this->m_g = g;

    this->m_n = n;
    if (L1) m_m1 = L1->dimensionOut();
    else m_m1 = n;
    if (L2) m_m2 = L2->dimensionOut();
    else m_m2 = n;
}

LinearOperator* FBProblem::getL1() {
    return m_L1;
}

LinearOperator* FBProblem::getL2() {
    return m_L2;
}

Matrix* FBProblem::getd1() {
    return m_d1;
}

Matrix* FBProblem::getd2() {
    return m_d2;
}

Function* FBProblem::getf1() {
    return m_f1;
}

Function* FBProblem::getf2() {
    return m_f2;
}

Function* FBProblem::getg() {
    return m_g;
}

Matrix* FBProblem::getlin() {
    return m_lin;
}

int FBProblem::getm1() {
    return m_m1;
}

int FBProblem::getm2() {
    return m_m2;
}

int FBProblem::getn() {
    return m_n;
}

FBProblem::~FBProblem() {
    // nothing to delete    
}
