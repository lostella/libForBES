/* 
 * File:   FunctionOntologyRegistry.h
 * Author: chung
 *
 * Created on October 28, 2015, 8:16 PM
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

#ifndef FUNCTIONONTOLOGYREGISTRY_H
#define	FUNCTIONONTOLOGYREGISTRY_H

#include "FunctionOntologicalClass.h"
#include "Function.h"
#include <iostream>

/**
 * 
 * \brief Registry of standard ontological classes.
 */
class FunctionOntologyRegistry {
public:

    static string nameSpace() {
        static string ns("fb");
        return ns;
    }

    static FunctionOntologicalClass function() {
        /*
         * NOTE:
         * The lifetime of function static variables begins the first time 
         * the program flow encounters the declaration and it ends at program 
         * termination.
         */
        static FunctionOntologicalClass generic_function("Function");
        return generic_function;
    }

    static FunctionOntologicalClass quadratic() {
        static FunctionOntologicalClass quad("Quadratic");
        quad.m_defines_f = true;
        quad.m_defines_conjugate = true;
        quad.m_defines_grad = true;
        quad.superClasses.push_back(function());
        return quad;
    }

    static FunctionOntologicalClass distance() {
        static FunctionOntologicalClass dist("Distance");
        dist.superClasses.push_back(function());
        return dist;
    }

    static FunctionOntologicalClass indicator() {
        static FunctionOntologicalClass ind("Indicator");
        ind.superClasses.push_back(function());
        return ind;
    }

    static FunctionOntologicalClass loss() {
        static FunctionOntologicalClass ind("LossFunction");
        ind.superClasses.push_back(function());
        return ind;
    }

    static FunctionOntologicalClass norm() {
        static FunctionOntologicalClass ind("Norm");
        ind.superClasses.push_back(function());
        return ind;
    }


private:
    FunctionOntologyRegistry();
    virtual ~FunctionOntologyRegistry();



};

#endif	/* FUNCTIONONTOLOGYREGISTRY_H */

