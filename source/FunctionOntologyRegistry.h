/* 
 * File:   FunctionOntologyRegistry.h
 * Author: Pantelis Sopasakis
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
#include <string>

/**
 * 
 * \brief Registry of standard ontological classes.
 */
class FunctionOntologyRegistry {
public:

    /**
     * Generic namespace of the libforbes function ontology
     * @return 
     */
    static std::string nameSpace();

    /**
     * Any proper convex closed extended real valued function \f$f:X\to\bar{\mathbb{R}}\f$
     * @return 
     */
    static FunctionOntologicalClass function();

    /**
     * Quadratic function, i.e., a function which can be written in the form
     * \f[
     *   f(x) = \frac{1}{2}x^{\top}Qx + q^{\top}x
     * \f]
     * @return 
     */
    static FunctionOntologicalClass quadratic();
    
    /**
     * Conjugate quadratic function, i.e., a function whose convex conjugate
     * is a quadratic function.
     * @return 
     */
    static FunctionOntologicalClass conj_quadratic();

    /**
     * A norm-distance from a nonempty convex closed set
     * @return 
     */
    static FunctionOntologicalClass distance();

    /**
     * Indicator of a nonempty convex closed set
     * @return 
     */
    static FunctionOntologicalClass indicator();

    /**
     * A loss function
     * @return 
     */
    static FunctionOntologicalClass loss();

    /**
     * A vector norm
     * @return 
     */
    static FunctionOntologicalClass norm();


private:


};

#endif	/* FUNCTIONONTOLOGYREGISTRY_H */

