/*
 * File:   TestOpAdjoint.h
 * Author: chung
 *
 * Created on Sep 15, 2015, 3:01:47 PM
 */

#ifndef TESTOPADJOINT_H
#define	TESTOPADJOINT_H
#define FORBES_TEST_UTILS

#include "ForBES.h"
#include <cppunit/extensions/HelperMacros.h>

class TestOpAdjoint : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestOpAdjoint);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallAdjoint);

    CPPUNIT_TEST_SUITE_END();

public:
    TestOpAdjoint();
    virtual ~TestOpAdjoint();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallAdjoint();

};

#endif	/* TESTOPADJOINT_H */

