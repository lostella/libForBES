/*
 * File:   TestFunctionOntologicalClass.h
 * Author: chung
 *
 * Created on Nov 2, 2015, 12:11:06 AM
 */

#ifndef TESTFUNCTIONONTOLOGICALCLASS_H
#define	TESTFUNCTIONONTOLOGICALCLASS_H

#define FORBES_TEST_UTILS
#include "ForBES.h"
#include <cppunit/extensions/HelperMacros.h>

class TestFunctionOntologicalClass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestFunctionOntologicalClass);

    CPPUNIT_TEST(testFunctionOntologicalClass);

    CPPUNIT_TEST_SUITE_END();

public:
    TestFunctionOntologicalClass();
    virtual ~TestFunctionOntologicalClass();
    void setUp();
    void tearDown();

private:
    void testFunctionOntologicalClass();

};

#endif	/* TESTFUNCTIONONTOLOGICALCLASS_H */

