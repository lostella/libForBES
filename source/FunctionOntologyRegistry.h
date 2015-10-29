/* 
 * File:   FunctionOntologyRegistry.h
 * Author: chung
 *
 * Created on October 28, 2015, 8:16 PM
 */

#ifndef FUNCTIONONTOLOGYREGISTRY_H
#define	FUNCTIONONTOLOGYREGISTRY_H

#include "FunctionOntologicalClass.h"
#include "Function.h"
#include <iostream>

class FunctionOntologyRegistry {
public:

    static string nameSpace() {
        static string ns("fb");
        return ns;
    }

    static FunctionOntologicalClass function() {
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
        static FunctionOntologicalClass ind("lossFunction");
        ind.superClasses.push_back(function());
        return ind;
    }
    
    static FunctionOntologicalClass norm() {
        static FunctionOntologicalClass ind("lossFunction");
        ind.superClasses.push_back(function());
        return ind;
    }


private:
    FunctionOntologyRegistry();
    virtual ~FunctionOntologyRegistry();



};

#endif	/* FUNCTIONONTOLOGYREGISTRY_H */

