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
 * 
 * \see FunctionOntologyRegistry
 */
class FunctionOntologicalClass {
public:

    /**
     * Creates a new ontological class given its name/unique identifier.
     * @param name Name of the ontological class.
     */
    FunctionOntologicalClass(string name);

    /**
     * Default destructor.
     */
    virtual ~FunctionOntologicalClass();

    /**
     * Unique name of the ontological class.
     * 
     * @return Unique name/identifier as a string
     */
    string getName() const;

    /** 
     * Whether this function type defines f(x).
     * @return <code>true</code> if f(x) is defined
     */
    bool defines_f() const {
        return m_defines_f;
    }

    /** 
     * Whether this function type defines the gradient of 
     * f at x, \f$\nabla f(x)\f$.
     * @return <code>true</code> if grad[f](x) is defined
     */
    bool defines_grad() const {
        return m_defines_grad;
    }

    /**
     * Whether this function type defines a conjugate \f$f^*(y)\f$
     * @return <code>true</code> if f*(y) is defined
     */
    bool defines_conjugate() const {
        return m_defines_conjugate;
    }

    /**
     * Whether this function type defines the gradient of its conjugate \f$\nabla f^*(y)\f$
     * @return <code>true</code> if grad[f*](y) is defined
     */
    bool defines_conjugate_grad() const {
        return m_defines_conjugate_grad;
    }

    /**
     * Whether this function defines a proximal \f$\mathrm{prox}_{\gamma f}(v)\f$
     * @return <code>true</code> if f*(x) is defined
     */
    bool defines_prox() const {
        return m_defines_prox;
    }

    list<FunctionOntologicalClass> getSuperclasses() const {
        return superClasses;
    }

    friend std::ostream& operator<<(std::ostream& os, const FunctionOntologicalClass& obj);

    void set_defines_conjugate(bool defines_conjugate) {
        m_defines_conjugate = defines_conjugate;
    }

    void set_defines_conjugate_grad(bool defines_conjugate_grad) {
        m_defines_conjugate_grad = defines_conjugate_grad;
    }

    void set_defines_f(bool defines_f) {
        m_defines_f = defines_f;
    }

    void set_defines_grad(bool defines_grad) {
        m_defines_grad = defines_grad;
    }

    void set_defines_prox(bool defines_prox) {
        m_defines_prox = defines_prox;
    }


private:

    FunctionOntologicalClass();

    friend class FunctionOntologyRegistry;

    string m_name;
    bool m_defines_f; /**< Whether this function type defines f(x).                                       */
    bool m_defines_grad; /**< Whether this function type defines the gradient of f at x, grad[f](x).         */
    bool m_defines_conjugate; /**< Whether this function type defines a conjugate f*(x)                           */
    bool m_defines_conjugate_grad; /**< Whether this function type defines the gradient of its conjugate grad[f*](x)   */
    bool m_defines_prox; /**< Whether this function defines a proximal prox_{gamma f}(v)                     */

    list<FunctionOntologicalClass> superClasses; /**< List of super-classes */

};



#endif	/* FUNCTIONONTOLOGICALCLASS_H */

