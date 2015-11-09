/* 
 * File:   FunctionOntologyRegistry.cpp
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

#include "FunctionOntologyRegistry.h"

FunctionOntologicalClass FunctionOntologyRegistry::distance() {
    static bool defines_conjugate = false;
    static bool defines_conjugate_grad = false;
    static bool defines_f = true;
    static bool defines_grad = false;
    static bool defines_prox = false;
    static bool defines_hess = false;
    static bool defines_hess_conj = false;
    static FunctionOntologicalClass dist(
            defines_conjugate,
            defines_conjugate_grad,
            defines_f,
            defines_grad,
            defines_prox,
            defines_hess,
            defines_hess_conj,
            "Distance",
            function());
    return dist;
}

FunctionOntologicalClass FunctionOntologyRegistry::function() {
    /*
     * NOTE:
     * The lifetime of function static variables begins the first time 
     * the program flow encounters the declaration and it ends at program 
     * termination.
     */
    static FunctionOntologicalClass generic_function("Function");
    return generic_function;
}

FunctionOntologicalClass FunctionOntologyRegistry::indicator() {
    static bool defines_conjugate = false;
    static bool defines_conjugate_grad = false;
    static bool defines_f = true;
    static bool defines_grad = false;
    static bool defines_prox = true;
    static bool defines_hess = false;
    static bool defines_hess_conj = false;
    static FunctionOntologicalClass dist(
            defines_conjugate,
            defines_conjugate_grad,
            defines_f,
            defines_grad,
            defines_prox,
            defines_hess,
            defines_hess_conj,
            "Indicator",
            function());
    return dist;
}

FunctionOntologicalClass FunctionOntologyRegistry::loss() {
    static bool defines_conjugate = false;
    static bool defines_conjugate_grad = false;
    static bool defines_f = true;
    static bool defines_grad = true;
    static bool defines_prox = false;
    static bool defines_hess = false;
    static bool defines_hess_conj = false;
    static FunctionOntologicalClass loss(
            defines_conjugate,
            defines_conjugate_grad,
            defines_f,
            defines_grad,
            defines_prox,
            defines_hess,
            defines_hess_conj,
            "LossFunction",
            function());
    return loss;
}

FunctionOntologicalClass FunctionOntologyRegistry::norm() {
    static bool defines_conjugate = true;
    static bool defines_conjugate_grad = false;
    static bool defines_f = true;
    static bool defines_grad = false;
    static bool defines_prox = true;
    static bool defines_hess = false;
    static bool defines_hess_conj = false;
    static FunctionOntologicalClass norm(
            defines_conjugate,
            defines_conjugate_grad,
            defines_f,
            defines_grad,
            defines_prox,
            defines_hess,
            defines_hess_conj,
            "Norm",
            function());
    return norm;
}

FunctionOntologicalClass FunctionOntologyRegistry::quadratic() {
    static bool defines_conjugate = true;
    static bool defines_conjugate_grad = true;
    static bool defines_f = true;
    static bool defines_grad = true;
    static bool defines_prox = false;
    static bool defines_hess = false;
    static bool defines_hess_conj = false;
    static FunctionOntologicalClass quad(
            defines_conjugate,
            defines_conjugate_grad,
            defines_f,
            defines_grad,
            defines_prox,
            defines_hess,
            defines_hess_conj,
            "Quadratic",
            function());
    return quad;
}

FunctionOntologicalClass FunctionOntologyRegistry::conj_quadratic() {
    static bool defines_conjugate = true;
    static bool defines_conjugate_grad = true;
    static bool defines_f = false;
    static bool defines_grad = false;
    static bool defines_prox = false;
    static bool defines_hess = false;
    static bool defines_hess_conj = false;
    static FunctionOntologicalClass quad(
            defines_conjugate,
            defines_conjugate_grad,
            defines_f,
            defines_grad,
            defines_prox,
            defines_hess,
            defines_hess_conj,
            "ConjugateQuadratic",
            function());
    return quad;
}

string FunctionOntologyRegistry::nameSpace() {
    static string ns("fb");
    return ns;
}

