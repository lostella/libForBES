/*
 * File:   TestFunctionOntologyRegistry.h
 * Author: chung
 *
 * Created on Nov 2, 2015, 12:21:37 AM
 */

#ifndef TESTFUNCTIONONTOLOGYREGISTRY_H
#define	TESTFUNCTIONONTOLOGYREGISTRY_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestFunctionOntologyRegistry : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestFunctionOntologyRegistry);

    CPPUNIT_TEST(testDistance);
    CPPUNIT_TEST(testFunction);
    CPPUNIT_TEST(testIndicator);
    CPPUNIT_TEST(testLoss);
    CPPUNIT_TEST(testNameSpace);
    CPPUNIT_TEST(testNorm);
    CPPUNIT_TEST(testQuadratic);

    CPPUNIT_TEST_SUITE_END();

public:
    TestFunctionOntologyRegistry();
    virtual ~TestFunctionOntologyRegistry();
    void setUp();
    void tearDown();

private:
    void testDistance();
    void testFunction();
    void testIndicator();
    void testLoss();
    void testNameSpace();
    void testNorm();
    void testQuadratic();

};

#endif	/* TESTFUNCTIONONTOLOGYREGISTRY_H */

