/*
 * File:   TestOntRegistry.h
 * Author: chung
 *
 * Created on Oct 29, 2015, 1:28:41 AM
 */

#ifndef TESTONTREGISTRY_H
#define	TESTONTREGISTRY_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestOntRegistry : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestOntRegistry);

    CPPUNIT_TEST(testOntologyRegistry);

    CPPUNIT_TEST_SUITE_END();

public:
    TestOntRegistry();
    virtual ~TestOntRegistry();
    void setUp();
    void tearDown();

private:
    void testOntologyRegistry();
};

#endif	/* TESTONTREGISTRY_H */

