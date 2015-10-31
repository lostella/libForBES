/* 
 * File:   FunctionOntologicalClass.cpp
 * Author: chung
 * 
 * Created on October 28, 2015, 8:07 PM
 */

#include "FunctionOntologicalClass.h"
#include "FunctionOntologyRegistry.h"

FunctionOntologicalClass::FunctionOntologicalClass() {
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

string FunctionOntologicalClass::getName() const {
    return m_name;
}

