/*
 * File:   TestSeparableSum.h
 * Author: chung
 *
 * Created on Nov 4, 2015, 12:10:44 AM
 */

#ifndef TESTSEPARABLESUM_H
#define	TESTSEPARABLESUM_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <vector>

#include <cppunit/extensions/HelperMacros.h>

class TestSeparableSum : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestSeparableSum);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallProx);

    CPPUNIT_TEST_SUITE_END();

public:
    TestSeparableSum();
    virtual ~TestSeparableSum();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallProx();
    
};

#endif	/* TESTSEPARABLESUM_H */

