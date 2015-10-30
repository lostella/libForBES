/*
 * File:   TestHuber.h
 * Author: chung
 *
 * Created on Oct 30, 2015, 2:19:35 AM
 */

#ifndef TESTHUBER_H
#define	TESTHUBER_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestHuber : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestHuber);

    CPPUNIT_TEST(testCall);

    CPPUNIT_TEST_SUITE_END();

public:
    TestHuber();
    virtual ~TestHuber();
    void setUp();
    void tearDown();

private:
    void testCall();

};

#endif	/* TESTHUBER_H */

