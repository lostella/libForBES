/*
 * File:   TestOntRegistry.cpp
 * Author: chung
 *
 * Created on Oct 29, 2015, 1:28:41 AM
 */

#include "TestOntRegistry.h"
#include "FunctionOntologicalClass.h"
#include "FunctionOntologyRegistry.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestOntRegistry);

TestOntRegistry::TestOntRegistry() {
}

TestOntRegistry::~TestOntRegistry() {
}

void TestOntRegistry::setUp() {
}

void TestOntRegistry::tearDown() {
}

void TestOntRegistry::testOntologyRegistry() {
    FunctionOntologicalClass myfun = FunctionOntologyRegistry::function();
    _ASSERT_NOT(myfun.defines_f());
    _ASSERT_NOT(myfun.defines_conjugate());
    _ASSERT_NOT(myfun.defines_grad());
    _ASSERT_NOT(myfun.defines_prox());
    _ASSERT_EQ(string("Function"), myfun.getName());
    list<FunctionOntologicalClass> superClasses = myfun.getSuperclasses();
    _ASSERT(superClasses.empty());
    
    myfun = FunctionOntologyRegistry::distance();
    superClasses = myfun.getSuperclasses();
    _ASSERT_NOT(superClasses.empty());
    
}

