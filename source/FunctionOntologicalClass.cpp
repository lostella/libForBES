/* 
 * File:   FunctionOntologicalClass.cpp
 * Author: Pantelis Sopasakis
 * 
 * Created on October 28, 2015, 8:07 PM
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

#include "FunctionOntologicalClass.h"
#include "FunctionOntologyRegistry.h"

FunctionOntologicalClass::FunctionOntologicalClass(
        bool does_define_conjugate,
        bool does_define_conjugate_grad,
        bool does_define_f,
        bool does_define_grad,
        bool does_define_prox,
        bool does_define_hessian,
        bool does_define_hessian_conj,
        string name,
        const FunctionOntologicalClass& super) :
m_name(name),
m_defines_f(does_define_f),
m_defines_grad(does_define_grad),
m_defines_conjugate(does_define_conjugate),
m_defines_conjugate_grad(does_define_conjugate_grad),
m_defines_prox(does_define_prox),
m_defines_hessian(does_define_hessian),
m_defines_hessian_conj(does_define_hessian_conj) {
    superClasses.push_back(super);
}

FunctionOntologicalClass::~FunctionOntologicalClass() {
}

FunctionOntologicalClass::FunctionOntologicalClass(string name) : m_name(name) {
    m_defines_f = false;
    m_defines_grad = false;
    m_defines_prox = false;
    m_defines_conjugate = false;
    m_defines_conjugate_grad = false;
}

//LCOV_EXCL_START

std::ostream& operator<<(std::ostream& os, const FunctionOntologicalClass& obj) {
    os << "Function class : " << obj.m_name << "\n";
    os << " * f()              : " << (obj.m_defines_f ? "YES" : "NO") << "\n";
    os << " * grad[f]()        : " << (obj.m_defines_grad ? "YES" : "NO") << "\n";
    os << " * hess[f]()        : " << (obj.m_defines_hessian ? "YES" : "NO") << "\n";
    os << " * f*()             : " << (obj.m_defines_conjugate ? "YES" : "NO") << "\n";
    os << " * grad[f*]()       : " << (obj.m_defines_conjugate_grad ? "YES" : "NO") << "\n";
    os << " * hess[f*]()       : " << (obj.m_defines_hessian_conj ? "YES" : "NO") << "\n";
    os << " * prox(gamma*f)()  : " << (obj.m_defines_prox ? "YES" : "NO") << "\n";
    os << "Super-classes... \n";
    list<FunctionOntologicalClass> li = obj.superClasses;
    size_t i = 1;
    for (list<FunctionOntologicalClass>::iterator it = li.begin(); it != li.end(); ++it) {
        FunctionOntologicalClass entry = *it;
        cout << " " << i << ". " << FunctionOntologyRegistry::nameSpace() << ":" << entry.getName() << "\n";
        i++;
    }
    return os;
}
//LCOV_EXCL_STOP

bool FunctionOntologicalClass::defines_conjugate() const {
    return m_defines_conjugate;
}

bool FunctionOntologicalClass::defines_conjugate_grad() const {
    return m_defines_conjugate_grad;
}

bool FunctionOntologicalClass::defines_f() const {
    return m_defines_f;
}

bool FunctionOntologicalClass::defines_grad() const {
    return m_defines_grad;
}

bool FunctionOntologicalClass::defines_prox() const {
    return m_defines_prox;
}

list<FunctionOntologicalClass> FunctionOntologicalClass::getSuperclasses() const {
    return superClasses;
}

void FunctionOntologicalClass::set_defines_conjugate(bool define_conjugate) {
    m_defines_conjugate = define_conjugate;
}

void FunctionOntologicalClass::set_defines_f(bool define_f) {
    m_defines_f = define_f;
}

void FunctionOntologicalClass::set_defines_conjugate_grad(bool define_conjugate_grad) {
    m_defines_conjugate_grad = define_conjugate_grad;
}

void FunctionOntologicalClass::set_defines_prox(bool define_prox) {
    m_defines_prox = define_prox;
}

void FunctionOntologicalClass::set_defines_grad(bool define_grad) {
    m_defines_grad = define_grad;
}

bool FunctionOntologicalClass::defines_hessian() const {
    return m_defines_hessian;
}

bool FunctionOntologicalClass::defines_hessian_conj() const {
    return m_defines_hessian_conj;
}

void FunctionOntologicalClass::set_defines_hessian(bool define_hessian) {
    m_defines_hessian = define_hessian;
}

void FunctionOntologicalClass::set_defines_hessian_conj(bool define_hessian_conj) {
    m_defines_hessian_conj = define_hessian_conj;
}

string FunctionOntologicalClass::getName() const {
    return this -> m_name;
}

