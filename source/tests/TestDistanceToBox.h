/*
 * File:   TestDistanceToBox.h
 * Author: chung
 *
 * Created on Oct 28, 2015, 7:28:27 PM
 */

#ifndef TESTDISTANCETOBOX_H
#define	TESTDISTANCETOBOX_H

#define FORBES_TEST_UTILS

#include "ForBES.h"
#include <cppunit/extensions/HelperMacros.h>

class TestDistanceToBox : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestDistanceToBox);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCall2);

    CPPUNIT_TEST_SUITE_END();

public:
    TestDistanceToBox();
    virtual ~TestDistanceToBox();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCall2();
    
};

#endif	/* TESTDISTANCETOBOX_H */

