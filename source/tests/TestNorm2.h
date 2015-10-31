/*
 * File:   TestNorm2.h
 * Author: chung
 *
 * Created on Oct 30, 2015, 9:53:13 PM
 */

#ifndef TESTNORM2_H
#define	TESTNORM2_H

#define FORBES_TEST_UTILS
#include "ForBES.h"

#include <cppunit/extensions/HelperMacros.h>

class TestNorm2 : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestNorm2);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testCallProx);
    CPPUNIT_TEST(testDualNorm);

    CPPUNIT_TEST_SUITE_END();

public:
    TestNorm2();
    virtual ~TestNorm2();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testCallProx();
    void testDualNorm();

};

#endif	/* TESTNORM2_H */

