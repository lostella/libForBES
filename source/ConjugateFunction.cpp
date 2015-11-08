/* 
 * File:   ConjugateFunction.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on November 7, 2015, 3:15 AM
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

#include "ConjugateFunction.h"

ConjugateFunction::ConjugateFunction(Function& funct) :
Function(), m_function(funct) {
    std::cout << "\n conjugate of " << funct.category().getName();
}

ConjugateFunction::~ConjugateFunction() {
}

int ConjugateFunction::call(Matrix& x, double& f) {
    return m_function.callConj(x, f);
}

int ConjugateFunction::call(Matrix& x, double& f, Matrix& grad) {
    return m_function.callConj(x, f, grad);
}

int ConjugateFunction::callConj(Matrix& x, double& f_star, Matrix& grad) {
    return m_function.call(x, f_star, grad);
}

int ConjugateFunction::callConj(Matrix& x, double& f_star) {
    return m_function.call(x, f_star);
}

int ConjugateFunction::callProx(Matrix& x, double gamma, Matrix& prox) {
    Matrix x_over_gamma = x;
    x_over_gamma *= (1.0/gamma);
    int status = m_function.callProx(x_over_gamma, 1.0/gamma, prox);
    prox *= gamma;
    for (size_t i = 0; i < x.getNrows(); i++) {
        prox.set(i, 0, x.get(i, 0) - prox.get(i, 0));
    }
    return status;
}

int ConjugateFunction::callProx(Matrix& x, double gamma, Matrix& prox, double& f_at_prox) {
    int status = callProx(x, gamma, prox);
    if (ForBESUtils::STATUS_OK != status){
        return status;
    }
    return call(prox, f_at_prox);
}

FunctionOntologicalClass ConjugateFunction::category() {
    FunctionOntologicalClass meta("Conjugate");
    FunctionOntologicalClass orig_meta = m_function.category();
    meta.set_defines_f(orig_meta.defines_conjugate());
    meta.set_defines_grad(orig_meta.defines_conjugate_grad());
    meta.set_defines_conjugate(orig_meta.defines_f());
    meta.set_defines_conjugate_grad(orig_meta.defines_grad());
    
    return meta;
}




