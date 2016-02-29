/*
 * File:   TestSumOfNorm2.h
 * Author: chung
 *
 * Created on Feb 29, 2016, 5:18:46 PM
 */

#ifndef TESTSUMOFNORM2_H
#define	TESTSUMOFNORM2_H

#include <cppunit/extensions/HelperMacros.h>

#define FORBES_TEST_UTILS
#include "ForBES.h"

class TestSumOfNorm2 : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TestSumOfNorm2);

    CPPUNIT_TEST(testCall);
    CPPUNIT_TEST(testDualNorm);
    CPPUNIT_TEST(testProx);

    CPPUNIT_TEST_SUITE_END();

public:
    TestSumOfNorm2();
    virtual ~TestSumOfNorm2();
    void setUp();
    void tearDown();

private:
    void testCall();
    void testDualNorm();
    void testProx();
    
};

#endif	/* TESTSUMOFNORM2_H */

