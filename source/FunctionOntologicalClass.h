/* 
 * File:   FunctionOntologicalClass.h
 * Author: chung
 *
 * Created on October 28, 2015, 8:07 PM
 */

#ifndef FUNCTIONONTOLOGICALCLASS_H
#define	FUNCTIONONTOLOGICALCLASS_H

#include <string>
#include <list>
#include <iostream>

using namespace std;

/**
 * Ontological class referring to a ForBES function.
 */
class FunctionOntologicalClass {
public:

    FunctionOntologicalClass();
    FunctionOntologicalClass(string name);
    virtual ~FunctionOntologicalClass();

    /**
     * Unique name of the ontological class.
     * 
     * @return Unique name/identifier as a string
     */
    string getName() const;    

    bool defines_conjugate() const {
        return m_defines_conjugate;
    }

    bool defines_f() const {
        return m_defines_f;
    }

    bool defines_grad() const {
        return m_defines_grad;
    }

    bool defines_prox() const {
        return m_defines_prox;
    }

    list<FunctionOntologicalClass> getSuperclasses() const {
        return superClasses;
    }
    
    friend std::ostream& operator<<(std::ostream& os, const FunctionOntologicalClass& obj);




private:
    friend class FunctionOntologyRegistry;

    string m_name;
    bool m_defines_f;
    bool m_defines_grad;
    bool m_defines_conjugate;
    bool m_defines_prox;

    list<FunctionOntologicalClass> superClasses;

};



#endif	/* FUNCTIONONTOLOGICALCLASS_H */

