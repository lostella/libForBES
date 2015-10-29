/* 
 * File:   FunctionOntologicalClass.h
 * Author: chung
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

#ifndef FUNCTIONONTOLOGICALCLASS_H
#define	FUNCTIONONTOLOGICALCLASS_H

#include <string>
#include <list>
#include <iostream>

using namespace std;

/**
 * \brief Ontological class referring to a ForBES function.
 */
class FunctionOntologicalClass {
    
public:
    
    /**
     * Creates a new ontological class given its name/unique identifier.
     * @param name Name of the ontological class.
     */
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
    
    FunctionOntologicalClass();
    
    friend class FunctionOntologyRegistry;

    string m_name;
    bool m_defines_f;
    bool m_defines_grad;
    bool m_defines_conjugate;
    bool m_defines_prox;

    list<FunctionOntologicalClass> superClasses;

};



#endif	/* FUNCTIONONTOLOGICALCLASS_H */

