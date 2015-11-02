/*
 * File:   TestFunctionOntologicalClass.cpp
 * Author: Pantelis Sopasakis
 *
 * Created on Nov 2, 2015, 12:11:07 AM
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

#include "TestFunctionOntologicalClass.h"
#include "FunctionOntologicalClass.h"
#include "ForBES.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestFunctionOntologicalClass);

TestFunctionOntologicalClass::TestFunctionOntologicalClass() {
}

TestFunctionOntologicalClass::~TestFunctionOntologicalClass() {
}

void TestFunctionOntologicalClass::setUp() {
}

void TestFunctionOntologicalClass::tearDown() {
}

void TestFunctionOntologicalClass::testFunctionOntologicalClass() {
    FunctionOntologicalClass ont("test");
    _ASSERT_EQ(string("test"), ont.getName());
    _ASSERT(ont.getSuperclasses().empty());
    
    _ASSERT_NOT(ont.defines_grad());
    ont.set_defines_grad(true);
    _ASSERT(ont.defines_grad());
    
    _ASSERT_NOT(ont.defines_f());
    ont.set_defines_f(true);
    _ASSERT(ont.defines_f());
    
    _ASSERT_NOT(ont.defines_conjugate());
    ont.set_defines_conjugate(true);
    _ASSERT(ont.defines_conjugate());
    
    _ASSERT_NOT(ont.defines_conjugate_grad());
    ont.set_defines_conjugate_grad(true);
    _ASSERT(ont.defines_conjugate_grad());
    
    _ASSERT_NOT(ont.defines_prox());
    ont.set_defines_prox(true);
    _ASSERT(ont.defines_prox());
    
}

