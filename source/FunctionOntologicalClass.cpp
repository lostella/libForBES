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

FunctionOntologicalClass::FunctionOntologicalClass(bool m_defines_conjugate, bool m_defines_conjugate_grad, bool m_defines_f, bool m_defines_grad, bool m_defines_prox, string m_name, const FunctionOntologicalClass& super) :
m_defines_conjugate(m_defines_conjugate),
m_defines_conjugate_grad(m_defines_conjugate_grad),
m_defines_f(m_defines_f),
m_defines_grad(m_defines_grad),
m_defines_prox(m_defines_prox),
m_name(m_name) {
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
    os << " * f*()             : " << (obj.m_defines_conjugate ? "YES" : "NO") << "\n";
    os << " * grad[f*]()       : " << (obj.m_defines_conjugate_grad ? "YES" : "NO") << "\n";
    os << " * prox(gamma*f)()  : " << (obj.m_defines_prox ? "YES" : "NO") << "\n";
    os << "Super-classes... \n";
    list<FunctionOntologicalClass> li = obj.superClasses;
    size_t i = 1;
    for (list<FunctionOntologicalClass>::iterator it = li.begin(); it != li.end(); it++) {
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

void FunctionOntologicalClass::set_defines_conjugate(bool defines_conjugate) {
    m_defines_conjugate = defines_conjugate;
}

void FunctionOntologicalClass::set_defines_f(bool defines_f) {
    m_defines_f = defines_f;
}

void FunctionOntologicalClass::set_defines_conjugate_grad(bool defines_conjugate_grad) {
    m_defines_conjugate_grad = defines_conjugate_grad;
}

void FunctionOntologicalClass::set_defines_prox(bool defines_prox) {
    m_defines_prox = defines_prox;
}

void FunctionOntologicalClass::set_defines_grad(bool defines_grad) {
    m_defines_grad = defines_grad;
}

string FunctionOntologicalClass::getName() const {
    return m_name;
}

