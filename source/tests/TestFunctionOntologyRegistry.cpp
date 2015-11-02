/*
 * File:   TestFunctionOntologyRegistry.cpp
 * Author: chung
 *
 * Created on Nov 2, 2015, 12:21:37 AM
 */

#include "TestFunctionOntologyRegistry.h"
#include "FunctionOntologyRegistry.h"


CPPUNIT_TEST_SUITE_REGISTRATION(TestFunctionOntologyRegistry);

TestFunctionOntologyRegistry::TestFunctionOntologyRegistry() {
}

TestFunctionOntologyRegistry::~TestFunctionOntologyRegistry() {
}

void TestFunctionOntologyRegistry::setUp() {
}

void TestFunctionOntologyRegistry::tearDown() {
}

void TestFunctionOntologyRegistry::testDistance() {
    FunctionOntologicalClass distance = FunctionOntologyRegistry::distance();
    _ASSERT(distance.defines_f());
    _ASSERT_NOT(distance.defines_conjugate());
    _ASSERT_EQ(string("Distance"), distance.getName());

}

void TestFunctionOntologyRegistry::testFunction() {
    FunctionOntologicalClass function = FunctionOntologyRegistry::function();
    _ASSERT_NOT(function.defines_f());
    _ASSERT_NOT(function.defines_grad());
    _ASSERT_NOT(function.defines_prox());
    _ASSERT_NOT(function.defines_conjugate());
    _ASSERT_EQ(string("Function"), function.getName());
}

void TestFunctionOntologyRegistry::testIndicator() {
    FunctionOntologicalClass indicator = FunctionOntologyRegistry::indicator();
    _ASSERT(indicator.defines_f());
    _ASSERT(indicator.defines_prox());
    _ASSERT_NOT(indicator.defines_conjugate());
    _ASSERT_NOT(indicator.defines_conjugate_grad());
    _ASSERT_NOT(indicator.defines_grad());
    _ASSERT_EQ(string("Indicator"), indicator.getName());
}

void TestFunctionOntologyRegistry::testLoss() {
    FunctionOntologicalClass loss = FunctionOntologyRegistry::loss();
    _ASSERT(loss.defines_f());
    _ASSERT_NOT(loss.defines_conjugate());
    _ASSERT_EQ(string("LossFunction"), loss.getName());
}

void TestFunctionOntologyRegistry::testNameSpace() {
    string result = FunctionOntologyRegistry::nameSpace();
    _ASSERT_EQ(string("fb"), result);
}

void TestFunctionOntologyRegistry::testNorm() {
    FunctionOntologicalClass norm = FunctionOntologyRegistry::norm();
    _ASSERT(norm.defines_f());
    _ASSERT(norm.defines_conjugate());
    _ASSERT(norm.defines_prox());
    _ASSERT_EQ(string("Norm"), norm.getName());
}

void TestFunctionOntologyRegistry::testQuadratic() {
    FunctionOntologicalClass quadratic = FunctionOntologyRegistry::quadratic();
    _ASSERT(quadratic.defines_f());
    _ASSERT(quadratic.defines_grad());
    _ASSERT(quadratic.defines_conjugate());
    _ASSERT(quadratic.defines_conjugate_grad());
    _ASSERT_EQ(string("Quadratic"), quadratic.getName());
}

