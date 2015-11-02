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
#include <iostream>

/**
 * 
 * \brief Registry of standard ontological classes.
 */
class FunctionOntologyRegistry {
public:

    static string nameSpace();

    static FunctionOntologicalClass function();

    static FunctionOntologicalClass quadratic();

    static FunctionOntologicalClass distance();

    static FunctionOntologicalClass indicator();

    static FunctionOntologicalClass loss();

    static FunctionOntologicalClass norm();


private:


};

#endif	/* FUNCTIONONTOLOGYREGISTRY_H */

