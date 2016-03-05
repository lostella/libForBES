/*
 * File:   TestFBStats.h
 * Author: chung
 *
 * Created on Mar 5, 2016, 4:53:48 PM
 */

#ifndef TESTFBSTATS_H
#define	TESTFBSTATS_H

#define FORBES_TEST_UTILS
#include "FBStats.h"
#include "ForBES.h"


#include <cppunit/extensions/HelperMacros.h>

class TestFBStats : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestFBStats);

    CPPUNIT_TEST(testDouble);
    CPPUNIT_TEST(testInteger);
    CPPUNIT_TEST(testVector);
    CPPUNIT_TEST(testNotExisting);
    

    CPPUNIT_TEST_SUITE_END();

public:
    TestFBStats();
    virtual ~TestFBStats();
    void setUp();
    void tearDown();

private:
    void testDouble();
    void testInteger();
    void testVector();
    void testNotExisting();
    
};

#endif	/* TESTFBSTATS_H */

